#!/usr/bin/env python3
import argparse
import csv
import os
import random
import sys
import multiprocessing as mp
from Bio.Seq import Seq

from pyPLP.thermodynamics import (
    calculate_tm,
    get_max_homopolymer,
    check_hairpin,
    compute_metrics
)
from pyPLP.multiplex.cross_dimer import calc_homodimer_dg, calc_heterodimer_dg
from pyPLP.plp_generation.segmented_backbone import resolve_restriction_site

def parse_args():
    parser = argparse.ArgumentParser(description="Generate thermodynamically sound Universal Detection Probes (Multiprocessed).")
    parser.add_argument("--length", type=int, default=18, help="Total length of the detection probe (default: 18)")
    parser.add_argument("--restriction-site", default="", help="Optional restriction enzyme cleavage site to embed (e.g., GCTGAGG)")
    parser.add_argument("--re-probe-offset", type=int, default=2, help="Padding on the 5' end of the embedded restriction site (default: 2)")
    parser.add_argument("--screen-against", nargs="*", default=[], help="Sequences or CSV files containing arms/PLPs to screen against to prevent cross-reactivity")
    parser.add_argument("--max-cross-dimer-dg", type=float, default=-5.0, help="Max allowed cross-reactivity dG (kcal/mol) with screen-against sequences (default: -5.0)")
    parser.add_argument("--num-candidates", type=int, default=5, help="Number of valid candidates to generate (default: 5)")
    parser.add_argument("--max-attempts", type=int, default=500000, help="Maximum number of generation attempts per core (default: 500000)")
    parser.add_argument("--threads", type=int, default=mp.cpu_count(), help="Number of CPU cores/threads to use (default: all)")
    return parser.parse_args()

def load_screen_sequences(inputs):
    seqs = set()
    for inp in inputs:
        if os.path.isfile(inp) and inp.lower().endswith(".csv"):
            with open(inp, "r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    if "Left_Arm" in row and row["Left_Arm"]:
                        seqs.add(row["Left_Arm"].strip().upper())
                    if "Right_Arm" in row and row["Right_Arm"]:
                        seqs.add(row["Right_Arm"].strip().upper())
                    if "PLP" in row and row["PLP"]:
                        seqs.add(row["PLP"].strip().upper())
        else:
            seqs.add(inp.strip().upper())
    return list(seqs)

def worker_search(worker_args):
    (
        length,
        re_seq,
        re_probe_offset,
        screen_seqs,
        max_cross_dimer_dg,
        max_attempts,
        worker_id
    ) = worker_args

    attempts = 0
    while attempts < max_attempts:
        attempts += 1
        
        if re_seq:
            right_pad_len = length - re_probe_offset - len(re_seq)
            if right_pad_len < 0:
                return {"error": f"Probe length {length} is too short to embed a {len(re_seq)}nt restriction site with a {re_probe_offset}nt offset."}
            left_pad = "".join(random.choices(["A", "C", "G", "T"], k=re_probe_offset))
            right_pad = "".join(random.choices(["A", "C", "G", "T"], k=right_pad_len))
            seq = left_pad + re_seq + right_pad
        else:
            seq = "".join(random.choices(["A", "C", "G", "T"], k=length))
            
        homo_len, homo_seq = get_max_homopolymer(seq)
        if homo_len >= 4:
            continue
            
        hp_found, hp_tm, _ = check_hairpin(seq)
        if hp_found == "Yes" and hp_tm > 30.0:
            continue
            
        tm = calculate_tm(seq)
        
        clash = False
        for s_seq in screen_seqs:
            dg = calc_heterodimer_dg(seq, s_seq)
            if dg < max_cross_dimer_dg:
                clash = True
                break
        if clash:
            continue
            
        self_dg = calc_homodimer_dg(seq)
        gc, _ = compute_metrics(Seq(seq))
        
        return {
            "seq": seq,
            "tm": tm,
            "gc": gc,
            "hp_tm": hp_tm,
            "self_dg": self_dg,
            "attempts": attempts
        }
    return None

def generate_probe():
    args = parse_args()
    
    re_seq = resolve_restriction_site(args.restriction_site) if args.restriction_site else ""
    screen_seqs = load_screen_sequences(args.screen_against)
    
    if screen_seqs:
        print(f"Loaded {len(screen_seqs)} sequence(s) to screen against cross-reactivity.")
        
    print(f"Generating up to {args.num_candidates} Universal Detection Probes utilizing {args.threads} CPU cores...")
    
    worker_args_list = []
    # We submit exactly as many tasks as candidates requested. Since they run in parallel, it'll scale across cores perfectly.
    for i in range(args.num_candidates):
        worker_args_list.append((
            args.length,
            re_seq,
            args.re_probe_offset,
            screen_seqs,
            args.max_cross_dimer_dg,
            args.max_attempts,
            i
        ))

    candidates = []
    total_attempts = 0
    with mp.Pool(processes=args.threads) as pool:
        results = pool.map(worker_search, worker_args_list)
        
        for res in results:
            if res is None:
                continue
            if "error" in res:
                print(f"Error: {res['error']}")
                sys.exit(1)
            total_attempts += res["attempts"]
            candidates.append(res)
        
    if not candidates:
        print(f"Failed to generate valid candidates after {args.max_attempts} attempts per core.")
        sys.exit(1)
        
    print("\n" + "="*80)
    print(f"{'Probe Sequence':<30} | {'Length':<6} | {'Tm (°C)':<7} | {'GC %':<6} | {'Hairpin':<7} | {'Self dG':<8}")
    print("-" * 80)
    for c in candidates:
        print(f"{c['seq']:<30} | {len(c['seq']):<6} | {c['tm']:<7.2f} | {c['gc']:<6.1f} | {c['hp_tm']:<7.2f} | {c['self_dg']:<8.2f}")
    print("="*80)
    print(f"Total cumulative attempts across cores to find {len(candidates)} candidates: {total_attempts}\n")

if __name__ == "__main__":
    generate_probe()
