#!/usr/bin/env python3
"""
design_plp.py: Standalone CLI and Python tool to design open, specific Padlock Probes (PLPs)
directly from custom Left Arm and Right Arm sequences outside the full pipeline.
"""

import os
import sys
import json
import csv
import argparse
import multiprocessing as mp
from pathlib import Path
from Bio.Seq import Seq

# Ensure pyPLP is importable
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

import pyPLP
from pyPLP.plp_generation.generator import generate_plp_for_arm
from pyPLP.plp_generation.worker import init_worker
from pyPLP.blast.discovery import discover_blast_dbs
from pyPLP.thermodynamics.metrics import compute_metrics
from pyPLP.thermodynamics.homopolymer import get_max_homopolymer
from pyPLP.thermodynamics.hairpin import check_hairpin


from pyPLP.plp_generation.structure_eval import eval_plp_structure
from pyPLP.plp_generation.vienna_plot import render_vienna_svg
from pyPLP.thermodynamics.mfold_detailed import eval_structure_mfold_detailed
from pyPLP.plp_generation.segmented_backbone import resolve_restriction_site
from pyPLP.crispr.designer import (
    determine_pam_and_rest_tail,
    determine_shared_pam_and_rest_tail,
    design_crispr_components,
    DEFAULT_CRRNA_REPEAT,
    CRISPR_REPORTERS,
)


def validate_arm_sequence(arm_seq: str, name: str) -> str:
    """Validates that an arm sequence consists of valid DNA characters."""
    clean = str(arm_seq).strip().upper()
    if not clean:
        raise ValueError(f"Error: {name} is empty.")
    invalid = set(clean) - set("ACGT")
    if invalid:
        raise ValueError(
            f"Error: {name} contains non-standard DNA nucleotide(s): {', '.join(sorted(invalid))}. "
            f"Only A, C, G, T are permitted."
        )
    return clean


def parse_arms_csv(csv_path: str) -> list:
    """
    Parses a CSV file containing multiple targets with Left_Arm and Right_Arm sequences.
    Returns a list of dicts: [{'target_id': str, 'left_arm': str, 'right_arm': str}, ...]
    """
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Input arms CSV '{csv_path}' not found.")

    targets = []
    with open(csv_path, "r", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        if not header:
            raise ValueError(f"CSV file '{csv_path}' is empty.")

        norm_header = [h.strip().lower() for h in header]
        left_col, right_col, id_col = None, None, None

        for idx, col in enumerate(norm_header):
            if col in ("left_arm", "leftarm", "left"):
                left_col = idx
            elif col in ("right_arm", "rightarm", "right"):
                right_col = idx
            elif col in ("", "target", "target_id", "target_name", "candidate_id", "id", "name"):
                id_col = idx

        if left_col is None or right_col is None:
            raise ValueError(
                f"CSV '{csv_path}' must contain 'Left_Arm' and 'Right_Arm' columns. Found columns: {header}"
            )

        for row_idx, row in enumerate(reader, start=1):
            if not row or not any(field.strip() for field in row):
                continue
            left_val = row[left_col].strip().upper() if left_col < len(row) else ""
            right_val = row[right_col].strip().upper() if right_col < len(row) else ""

            if not left_val or not right_val:
                continue

            clean_left = validate_arm_sequence(left_val, f"Row {row_idx} Left Arm")
            clean_right = validate_arm_sequence(right_val, f"Row {row_idx} Right Arm")

            target_id = ""
            if id_col is not None and id_col < len(row):
                target_id = row[id_col].strip()
            if not target_id:
                target_id = f"Target_{len(targets) + 1}"

            targets.append({
                "target_id": target_id,
                "left_arm": clean_left,
                "right_arm": clean_right,
            })

    if not targets:
        raise ValueError(f"No valid target arm pairs found in '{csv_path}'.")

    return targets


def design_custom_plp(
    left_arm: str,
    right_arm: str,
    job_name: str,
    count: int = None,
    backbone_len: int = None,
    backbone_length_factor: int = None,
    output_csv: str = None,
    structures_dir: str = None,
    blastdb_dir: str = "blastdb",
    skip_blast: bool = False,
    threads: int = None,
    config_file: str = "config.json",
    keep_ct_ps: bool = False,
    probe_prefix: str = "Custom_PLP",
    universal_rca_primer: str = None,
    universal_probe: str = None,
    spacer1_len: int = 4,
    spacer2_len: int = 4,
    rca_primer_binding_len: int = 20,
    probe_binding_len: int = 20,
    spacer3_len: int = 0,
    rest_len: int = 4,
    restriction_site: str = None,
    screen_against: list = None,
    max_cross_dimer_dg: float = None,
    crispr_mode: bool = False,
    crrna_repeat: str = DEFAULT_CRRNA_REPEAT,
    spacer_len: int = 20,
    cis_probe_len: int = 30,
    generate_cis_probe: bool = False,
) -> list:
    """
    Designs confirmed open, non-cross-reactive padlock probes for a given Left and Right arm pair,
    with segmented LFA backbone architecture and automated RCA Primer & FAM Probe generation.
    All outputs (probes CSV, LFA order sheet, and structure plots) are saved inside results/<job_name>/.
    """
    clean_left = validate_arm_sequence(left_arm, "Left Arm")
    clean_right = validate_arm_sequence(right_arm, "Right Arm")

    if not job_name or not str(job_name).strip():
        raise ValueError("Error: job_name (-j / --job-name) is required.")

    raw_job = str(job_name).strip()
    job_path = Path(raw_job)
    if not job_path.is_absolute():
        parts = job_path.parts
        if not parts or parts[0] != "results":
            job_path = Path("results") / job_path
    job_dir = str(job_path.resolve())
    os.makedirs(job_dir, exist_ok=True)

    # Load base configuration if available
    cfg = {}
    if config_file and os.path.exists(config_file):
        try:
            with open(config_file, "r") as f:
                cfg = json.load(f)
        except Exception as e:
            print(f"[Notice] Could not load '{config_file}': {e}. Using defaults.")

    arm_cfg = cfg.get("arm_design", {})
    plp_cfg = cfg.get("plp_generation", {})

    target_count = int(count if count is not None else plp_cfg.get("plp_per_arm", 5))

    # Determine CSV output path
    if not output_csv:
        final_csv_path = os.path.join(job_dir, "designed_plp_probes.csv")
    elif os.path.isabs(output_csv) or os.path.dirname(output_csv):
        final_csv_path = os.path.abspath(output_csv)
    else:
        final_csv_path = os.path.join(job_dir, output_csv)

    # Determine structures directory path
    configured_struct_dir = structures_dir or plp_cfg.get("structures_dir") or "plp_structures"
    if os.path.isabs(configured_struct_dir) or os.path.dirname(configured_struct_dir):
        final_struct_dir = os.path.abspath(configured_struct_dir)
    else:
        final_struct_dir = os.path.join(job_dir, configured_struct_dir)

    os.makedirs(final_struct_dir, exist_ok=True)

    len_l = len(clean_left)
    len_r = len(clean_right)
    avg_arm = (len_l + len_r) // 2

    # Resolve settings with CLI / explicit overrides
    num_threads = int(threads or plp_cfg.get("threads") or arm_cfg.get("threads") or os.cpu_count() or 32)

    arch_cfg = plp_cfg.get("backbone_architecture", {})
    s1_len = int(spacer1_len if spacer1_len is not None else arch_cfg.get("spacer1_len", 4))
    rca_len = int(rca_primer_binding_len if rca_primer_binding_len is not None else arch_cfg.get("rca_primer_binding_len", 20))
    s2_len = int(spacer2_len if spacer2_len is not None else arch_cfg.get("spacer2_len", 4))
    pr_len = int(probe_binding_len if probe_binding_len is not None else arch_cfg.get("probe_binding_len", 20))
    s3_len = int(spacer3_len if spacer3_len is not None else arch_cfg.get("spacer3_len", 0))
    r_len = int(rest_len if rest_len is not None else arch_cfg.get("rest_len", 4))
    univ_rca = (universal_rca_primer if universal_rca_primer is not None else arch_cfg.get("universal_rca_primer_seq", "")).strip()
    if univ_rca:
        rca_len = len(univ_rca)
        
    univ_probe = (universal_probe if universal_probe is not None else arch_cfg.get("universal_probe_seq", "")).strip()
    if univ_probe:
        pr_len = len(univ_probe)
        
    min_tm_filter = float(arm_cfg.get("min_tm", 40.0))

    re_site_str = str(restriction_site if restriction_site is not None else arch_cfg.get("restriction_site", "")).strip()
    from pyPLP.plp_generation.segmented_backbone import resolve_restriction_site
    re_seq = resolve_restriction_site(re_site_str) if re_site_str else ""
    re_len = len(re_seq)
    
    re_mode = arch_cfg.get("re_mode", "standalone")
    re_probe_offset = int(arch_cfg.get("re_probe_offset", 2))
    actual_re_len = re_len if re_mode == "standalone" else 0

    pam_info = {}
    if crispr_mode:
        pam_info = determine_pam_and_rest_tail(clean_right)
        rest_tail = pam_info["rest_tail"]
        if len(rest_tail) >= r_len:
            r_len = len(rest_tail) + 1

    arch_bb_len = s1_len + rca_len + s2_len + pr_len + s3_len + actual_re_len + r_len

    eff_bb_factor = backbone_length_factor if backbone_length_factor is not None else plp_cfg.get("backbone_length_factor")
    if backbone_len is not None:
        target_bb_len = int(backbone_len)
    elif plp_cfg.get("backbone_length"):
        target_bb_len = int(plp_cfg.get("backbone_length"))
    elif eff_bb_factor is not None and str(eff_bb_factor).strip():
        target_bb_len = int(eff_bb_factor) * avg_arm
    else:
        target_bb_len = arch_bb_len

    if target_bb_len < arch_bb_len:
        actual_bb_len = arch_bb_len
    else:
        actual_bb_len = target_bb_len
        r_len += (actual_bb_len - arch_bb_len)
        arch_bb_len = actual_bb_len

    # Ensure r_len is at least the length of the engineered rest_tail in CRISPR mode
    if crispr_mode and pam_info.get("rest_tail"):
        req_tail_len = len(pam_info["rest_tail"])
        if r_len < req_tail_len:
            actual_bb_len += (req_tail_len - r_len)
            r_len = req_tail_len
            arch_bb_len = actual_bb_len

    import csv
    parsed_screen_seqs = set()
    if screen_against:
        for inp in screen_against:
            if os.path.isfile(inp) and inp.lower().endswith(".csv"):
                with open(inp, "r") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        for col in ["Left_Arm", "Right_Arm", "PLP"]:
                            if col in row and row[col]:
                                parsed_screen_seqs.add(row[col].strip().upper())
            else:
                parsed_screen_seqs.add(inp.strip().upper())

    # Do not screen the PLP against its own Left and Right arms
    parsed_screen_seqs.discard(clean_left)
    parsed_screen_seqs.discard(clean_right)

    # Inherit or fallback thermodynamic thresholds
    generation_cfg = {
        "screen_against": list(parsed_screen_seqs),
        "max_cross_dimer_dg": float(max_cross_dimer_dg) if max_cross_dimer_dg is not None else float(plp_cfg.get("max_cross_dimer_dg", -5.0)),
        "plp_per_arm": target_count,
        "threads": num_threads,
        "backbone_length": actual_bb_len,
        "maxHomopolymer": max(int(plp_cfg.get("maxHomopolymer", arm_cfg.get("maxHomopolymer", 4))), 4 if crispr_mode else 1),
        "max_hairpin_tm": float(plp_cfg.get("max_hairpin_tm", arm_cfg.get("max_hairpin_tm", 30.0))),
        "min_plp_mfe_dg": float(plp_cfg.get("min_plp_mfe_dg", -6.0)),
        "max_paired_arm_bases": int(plp_cfg.get("max_paired_arm_bases", 4)),
        "max_paired_total_bases": int(plp_cfg.get("max_paired_total_bases", 16)),
        "blastMinIdentity": float(plp_cfg.get("blastMinIdentity", arm_cfg.get("blastMinIdentity", 70.0))),
        "maxAllowedHits": plp_cfg.get("maxAllowedHits", arm_cfg.get("maxAllowedHits", 0)),
        "keep_ct_ps": keep_ct_ps or bool(plp_cfg.get("keep_ct_ps", False)),
        "backbone_architecture": {
            "spacer1_len": s1_len,
            "rca_primer_binding_len": rca_len,
            "spacer2_len": s2_len,
            "probe_binding_len": pr_len,
            "spacer3_len": s3_len,
            "rest_len": r_len,
            "restriction_site": re_site_str,
            "re_mode": re_mode,
            "re_probe_offset": re_probe_offset,
            "universal_rca_primer_seq": univ_rca,
            "universal_probe_seq": univ_probe,
            "min_rca_primer_tm": min_tm_filter,
            "min_probe_tm": min_tm_filter,
            "crispr_mode": crispr_mode,
            "rest_tail": pam_info.get("rest_tail", "") if crispr_mode else "",
        },
    }

    # Discover BLAST background databases
    blast_dbs = {}
    if not skip_blast and blastdb_dir and os.path.exists(blastdb_dir):
        blast_dbs = discover_blast_dbs(blastdb_dir)

    # Compute arm biophysical profiles
    l_gc, l_tm = compute_metrics(Seq(clean_left))
    r_gc, r_tm = compute_metrics(Seq(clean_right))
    l_homo_len, l_homo_seq = get_max_homopolymer(clean_left)
    r_homo_len, r_homo_seq = get_max_homopolymer(clean_right)
    l_hp_found, l_hp_tm, _ = check_hairpin(clean_left)
    r_hp_found, r_hp_tm, _ = check_hairpin(clean_right)

    print("=" * 80)
    print("STANDALONE PADLOCK PROBE (PLP) DESIGNER" + (" [CRISPR Cas12a Mode]" if crispr_mode else ""))
    print("=" * 80)
    print(f"Job Name (Output Dir)  : {job_dir}")
    print(f"Left Arm  (5'->3')     : {clean_left} ({len_l} nt, GC: {l_gc}%, Tm: {l_tm}°C, Hairpin: {l_hp_found} [{l_hp_tm}°C])")
    print(f"Right Arm (5'->3')     : {clean_right} ({len_r} nt, GC: {r_gc}%, Tm: {r_tm}°C, Hairpin: {r_hp_found} [{r_hp_tm}°C])")
    if crispr_mode:
        print(f"CRISPR Cas Mode        : ENABLED")
        print(f"Target PAM Motif       : {pam_info.get('pam_seq', '[A/C/G]TT[A/C/G]')} ([A/C/G]-T-T-[A/C/G] / {pam_info.get('pam_location', '')})")
        print(f"Engineered Rest Tail   : '{pam_info.get('rest_tail', '')}' (appended at Rest 3'-end)")
        print(f"Cis-Probe Generation   : {'Enabled' if generate_cis_probe else 'Disabled (Default)'}")
        print(f"crRNA Repeat Sequence  : {crrna_repeat} ({len(crrna_repeat)} nt)")
        print(f"crRNA Spacer Length    : {spacer_len} nt")
        print(f"CRISPR FQ Reporter     : 5'-FAM-TTATTATT-BHQ1-3'")
        print(f"CRISPR FB Reporter     : 5'-FAM-TTATTATT-Biotin-3'")
    print(f"Backbone Architecture  : Spacer1({s1_len}nt) + RCA({rca_len}nt) + Spacer2({s2_len}nt) + Probe({pr_len}nt) + Spacer3({s3_len}nt) + RE({actual_re_len}nt) + Rest({r_len}nt) = {actual_bb_len} nt")
    print(f"Total Intact PLP Length: {len_l + actual_bb_len + len_r} nt")
    print(f"Target Probes Count    : {target_count}")
    print(f"Openness Criteria      : ViennaRNA dG >= {generation_cfg['min_plp_mfe_dg']} kcal/mol, Arm Paired <= {generation_cfg['max_paired_arm_bases']} nt, Total Paired <= {generation_cfg['max_paired_total_bases']} nt")
    print(f"Minimum Tm Criteria    : RCA Primer Tm >= {generation_cfg['backbone_architecture']['min_rca_primer_tm']}°C, Detection Probe Tm >= {generation_cfg['backbone_architecture']['min_probe_tm']}°C")
    print(f"Active CPU Cores       : {num_threads}")
    print(f"Probes CSV Destination : {final_csv_path}")
    print(f"Structures Directory   : {final_struct_dir}")
    print(f"BLAST Background DBs   : {', '.join(blast_dbs.keys()) if blast_dbs else ('Skipped (--skip-blast)' if skip_blast else 'None found')}")
    print("-" * 80 + "\n")

    total_counter = mp.Value("L", 0)
    lock = mp.Lock()
    stop_event = mp.Event()

    with mp.Pool(
        processes=num_threads,
        initializer=init_worker,
        initargs=(total_counter, lock, stop_event),
    ) as pool:
        records = generate_plp_for_arm(
            left_arm=clean_left,
            right_arm=clean_right,
            arm_meta={"Candidate_ID": probe_prefix},
            cfg=generation_cfg,
            structures_dir=final_struct_dir,
            blast_dbs=blast_dbs,
            pool=pool,
            total_counter=total_counter,
            lock=lock,
            stop_event=stop_event,
        )

    for idx, r in enumerate(records, start=1):
        r["plp_id"] = f"{probe_prefix}_{idx}"

    if crispr_mode:
        for rec in records:
            rec["crispr"] = design_crispr_components(
                left_arm=clean_left,
                right_arm=clean_right,
                backbone=rec["bb"],
                crrna_repeat=crrna_repeat,
                spacer_len=spacer_len,
                cis_probe_len=cis_probe_len,
                generate_cis_probe=generate_cis_probe,
            )

    headers = [
        "PLP_ID",
        "Full_PLP_Sequence",
        "Total_Length",
        "Full_PLP_GC%",
        "Full_PLP_Tm",
        "Left_Arm",
        "Backbone",
        "Right_Arm",
        "Backbone_Length",
        "Backbone_GC%",
        "Backbone_Tm",
        "Backbone_Max_Homopolymer",
        "Backbone_Hairpin_Found",
        "Backbone_Hairpin_Tm",
        "Spacer1_Seq",
        "RCA_Primer_Binding_Region",
        "Spacer2_Seq",
        "Probe_Binding_Region",
        "Rest_Backbone_Seq",
        "RCA_Primer_Seq_5to3",
        "RCA_Primer_5_Mod",
        "RCA_Primer_Tm",
        "RCA_Primer_GC%",
        "Detection_Probe_Seq_5to3",
        "Detection_Probe_5_Mod",
        "Detection_Probe_Tm",
        "Detection_Probe_GC%",
        "Vienna_Full_PLP_Structure",
        "Vienna_Full_PLP_dG_kcal_mol",
        "mfold_Full_PLP_Structure",
        "mfold_Full_PLP_dG_kcal_mol",
        "Left_Arm_Paired_Bases",
        "Right_Arm_Paired_Bases",
        "Total_Paired_Bases",
        "Vienna_Structure_Plot",
        "mfold_Structure_Plot",
    ]
    if crispr_mode:
        headers.extend([
            "CRISPR_PAM_Motif",
            "CRISPR_PAM_Location",
        ])
        if generate_cis_probe:
            headers.extend([
                "CRISPR_Cis_Probe_Seq_5to3",
                "CRISPR_Cis_Probe_Length",
                "CRISPR_Cis_Probe_Tm",
                "CRISPR_Cis_Probe_GC%",
            ])
        headers.extend([
            "CRISPR_crRNA_Full_Seq_5to3",
            "CRISPR_crRNA_Repeat",
            "CRISPR_crRNA_Spacer",
            "CRISPR_crRNA_Spacer_Length",
        ])
    for db_name in blast_dbs:
        headers.append(f"Backbone_Hits_{db_name}")

    os.makedirs(os.path.dirname(os.path.abspath(final_csv_path)) or ".", exist_ok=True)
    with open(final_csv_path, "w", newline="") as f_out:
        writer = csv.writer(f_out)
        writer.writerow(headers)
        for rec in records:
            row = [
                rec["plp_id"],
                rec["full_plp"],
                len(rec["full_plp"]),
                rec["full_gc"],
                rec["full_tm"],
                clean_left,
                rec["bb"],
                clean_right,
                len(rec["bb"]),
                rec["bb_gc"],
                rec["bb_tm"],
                rec["bb_homo"],
                rec["bb_hp_found"],
                rec["bb_hp_tm"],
                rec.get("Spacer1_Seq", ""),
                rec.get("RCA_Primer_Binding_Region", ""),
                rec.get("Spacer2_Seq", ""),
                rec.get("Probe_Binding_Region", ""),
                rec.get("Rest_Backbone_Seq", ""),
                rec.get("RCA_Primer_Seq_5to3", ""),
                rec.get("RCA_Primer_5_Mod", "Biotin"),
                rec.get("RCA_Primer_Tm", ""),
                rec.get("RCA_Primer_GC%", ""),
                rec.get("Detection_Probe_Seq_5to3", ""),
                rec.get("Detection_Probe_5_Mod", "6-FAM"),
                rec.get("Detection_Probe_Tm", ""),
                rec.get("Detection_Probe_GC%", ""),
                rec["dot_bracket"],
                rec["mfe_dg"],
                rec["mfold_db"],
                rec["mfold_dg"],
                rec["paired_left"],
                rec["paired_right"],
                rec["paired_total"],
                rec["vienna_plot"],
                rec["mfold_plot"],
            ]
            if crispr_mode:
                c = rec.get("crispr", {})
                row.extend([
                    c.get("pam_seq", ""),
                    c.get("pam_location", ""),
                ])
                if generate_cis_probe:
                    row.extend([
                        c.get("cis_probe_seq", ""),
                        c.get("cis_probe_len", ""),
                        c.get("cis_probe_tm", ""),
                        c.get("cis_probe_gc", ""),
                    ])
                row.extend([
                    c.get("crrna_full_seq", ""),
                    c.get("crrna_repeat", ""),
                    c.get("spacer_seq_rna", ""),
                    c.get("spacer_len", ""),
                ])
            row.extend(rec.get("blast_hits", []))
            writer.writerow(row)

    # Generate dedicated LFA order sheet CSV
    lfa_csv = os.path.join(job_dir, "lfa_rca_oligos.csv")
    lfa_headers = [
        "Oligo_Name",
        "Sequence_5to3",
        "5_Modification",
        "Length_nt",
        "GC%",
        "Tm_C",
        "Role",
        "Associated_PLP_ID",
    ]
    with open(lfa_csv, "w", newline="") as f_lfa:
        lfa_writer = csv.writer(f_lfa)
        lfa_writer.writerow(lfa_headers)
        for rec in records:
            pid = rec["plp_id"]
            lfa_writer.writerow([
                pid,
                rec["full_plp"],
                "5'-Phos",
                len(rec["full_plp"]),
                rec["full_gc"],
                rec["full_tm"],
                "Padlock Probe (PLP)",
                pid,
            ])
            if rec.get("RCA_Primer_Seq_5to3"):
                lfa_writer.writerow([
                    f"{pid}_RCA_Primer",
                    rec["RCA_Primer_Seq_5to3"],
                    "5'-Biotin",
                    len(rec["RCA_Primer_Seq_5to3"]),
                    rec.get("RCA_Primer_GC%", ""),
                    rec.get("RCA_Primer_Tm", ""),
                    "RCA Amplification Primer",
                    pid,
                ])
            if rec.get("Detection_Probe_Seq_5to3"):
                lfa_writer.writerow([
                    f"{pid}_Detection_Probe",
                    rec["Detection_Probe_Seq_5to3"],
                    "5'-6-FAM",
                    len(rec["Detection_Probe_Seq_5to3"]),
                    rec.get("Detection_Probe_GC%", ""),
                    rec.get("Detection_Probe_Tm", ""),
                    "LFA Fluorescent Detection Probe",
                    pid,
                ])
        if crispr_mode:
            for rep in CRISPR_REPORTERS:
                lfa_writer.writerow([
                    rep["oligo_name"],
                    rep["sequence_5to3"],
                    rep["modifications"],
                    rep["length_nt"],
                    rep["gc_percent"],
                    rep["tm_c"],
                    rep["role"],
                    "Universal",
                ])

    # Generate dedicated CRISPR Cas12a order sheet CSV if crispr_mode
    crispr_csv = None
    if crispr_mode:
        crispr_csv = os.path.join(job_dir, "crispr_oligos.csv")
        crispr_headers = [
            "Oligo_Name",
            "Sequence_5to3",
            "Modifications",
            "Type",
            "Length_nt",
            "GC%",
            "Tm_C",
            "Role",
            "Associated_PLP_ID",
        ]
        with open(crispr_csv, "w", newline="") as f_cr:
            cr_writer = csv.writer(f_cr)
            cr_writer.writerow(crispr_headers)
            for rec in records:
                pid = rec["plp_id"]
                c = rec.get("crispr", {})
                cr_writer.writerow([
                    pid,
                    rec["full_plp"],
                    "5'-Phosphorylated",
                    "DNA",
                    len(rec["full_plp"]),
                    rec["full_gc"],
                    rec["full_tm"],
                    "Padlock Probe (PLP, 5' Phosphorylated for ligation)",
                    pid,
                ])
                if c.get("cis_probe_seq"):
                    cr_writer.writerow([
                        f"{pid}_Cis_Probe",
                        c.get("cis_probe_seq", ""),
                        "None",
                        "DNA",
                        c.get("cis_probe_len", ""),
                        c.get("cis_probe_gc", ""),
                        c.get("cis_probe_tm", ""),
                        f"CRISPR Cas12a Cis-Probe Activator ({c.get('pam_location', '')})",
                        pid,
                    ])
                cr_writer.writerow([
                    f"{pid}_crRNA",
                    c.get("crrna_full_seq", ""),
                    "None (RNA)",
                    "RNA",
                    c.get("crrna_full_len", len(c.get("crrna_full_seq", ""))),
                    c.get("spacer_gc", ""),
                    c.get("spacer_tm", ""),
                    f"crRNA Guide (Repeat: {len(c.get('crrna_repeat', ''))}nt, Spacer: {c.get('spacer_len', '')}nt)",
                    pid,
                ])
            # Always append Universal CRISPR Cas12a Reporter Probes
            for rep in CRISPR_REPORTERS:
                cr_writer.writerow([
                    rep["oligo_name"],
                    rep["sequence_5to3"],
                    rep["modifications"],
                    rep["type"],
                    rep["length_nt"],
                    rep["gc_percent"],
                    rep["tm_c"],
                    rep["role"],
                    "Universal",
                ])

    print("\n" + "=" * 80)
    print(f"SUCCESS: Generated {len(records)} confirmed open padlock probes!")
    print(f"Job Output Directory  -> {job_dir}")
    print(f"Exported probes CSV   -> {final_csv_path}")
    print(f"Exported LFA oligos   -> {lfa_csv}")
    if crispr_csv:
        print(f"Exported CRISPR oligos-> {crispr_csv}")
    print(f"Vienna SVG plots      -> {os.path.abspath(os.path.join(final_struct_dir, 'vienna'))}/")
    print(f"DNA mfold PNG plots   -> {os.path.abspath(os.path.join(final_struct_dir, 'mfold'))}/")
    print("=" * 80 + "\n")

    # Display clean terminal summary table
    print(f"{'PLP ID':<16} {'Length':<8} {'Vienna dG':<11} {'mfold dG':<11} {'RCA Primer (Tm)':<22} {'FAM Probe (Tm)':<20}")
    print("-" * 92)
    for rec in records:
        rca_info = f"{rec.get('RCA_Primer_Seq_5to3', '')[:8]}.. ({rec.get('RCA_Primer_Tm', '')}°C)"
        probe_info = f"{rec.get('Detection_Probe_Seq_5to3', '')[:8]}.. ({rec.get('Detection_Probe_Tm', '')}°C)"
        print(
            f"{rec['plp_id']:<16} "
            f"{len(rec['full_plp']):<8} "
            f"{rec['mfe_dg']:<11.2f} "
            f"{rec['mfold_dg']:<11.2f} "
            f"{rca_info:<22} "
            f"{probe_info:<20}"
        )
    print("-" * 92 + "\n")

    if crispr_mode:
        print("CRISPR CAS12A ASSAY OLIGOS & crRNA GUIDES:")
        if generate_cis_probe:
            print(f"{'PLP ID':<16} {'PAM':<6} {'Cis Probe (5\'->3\')':<34} {'crRNA Spacer (5\'->3\')':<24}")
            print("-" * 84)
            for rec in records:
                c = rec.get("crispr", {})
                print(f"{rec['plp_id']:<16} {c.get('pam_seq', ''):<6} {c.get('cis_probe_seq', ''):<34} {c.get('spacer_seq_rna', ''):<24}")
            print("-" * 84)
        else:
            print(f"{'PLP ID':<16} {'PAM':<6} {'crRNA Spacer (5\'->3\')':<24} {'PAM Location':<30}")
            print("-" * 80)
            for rec in records:
                c = rec.get("crispr", {})
                print(f"{rec['plp_id']:<16} {c.get('pam_seq', ''):<6} {c.get('spacer_seq_rna', ''):<24} {c.get('pam_location', ''):<30}")
            print("-" * 80)
        print("Universal CRISPR Cas12a Reporter Probes:")
        for rep in CRISPR_REPORTERS:
            print(f"  * {rep['oligo_name']:<26}: {rep['sequence_5to3']} ({rep['role']})")
        print("-" * 80 + "\n")

    return records


def design_shared_backbone_from_csv(
    arms_csv: str,
    job_name: str,
    count: int = 100,
    backbone_len: int = None,
    backbone_length_factor: int = None,
    output_csv: str = None,
    structures_dir: str = None,
    blastdb_dir: str = "blastdb",
    skip_blast: bool = False,
    threads: int = None,
    config_file: str = "config.json",
    keep_ct_ps: bool = False,
    probe_prefix: str = "Universal_PLP",
    universal_rca_primer: str = None,
    universal_probe: str = None,
    spacer1_len: int = 4,
    spacer2_len: int = 4,
    rca_primer_binding_len: int = 20,
    probe_binding_len: int = 20,
    spacer3_len: int = 0,
    rest_len: int = 4,
    restriction_site: str = None,
    screen_against: list = None,
    max_cross_dimer_dg: float = None,
    crispr_mode: bool = False,
    crrna_repeat: str = "UAAUUUCUACUAAGUGUAGAU",
    spacer_len: int = 20,
    cis_probe_len: int = 28,
    generate_cis_probe: bool = False,
) -> list:
    """
    Designs a shared/universal backbone across multiple target arm pairs provided in a CSV file.
    1. Generates candidate backbones on Target_1 (seed/reference).
    2. Grafts each backbone onto Target_2 through Target_N.
    3. Computes secondary structure folding dG and arm base pairing for all targets.
    4. Ranks and sorts the backbones by cross-target openness (all_pass, pass_count, worst_dg, avg_dg).
    5. Exports ranked candidate backbones, full probes for all targets, and LFA order sheets.
    """
    targets = parse_arms_csv(arms_csv)
    num_targets = len(targets)

    if not job_name or not str(job_name).strip():
        raise ValueError("Error: job_name (-j / --job-name) is required.")

    raw_job = str(job_name).strip()
    job_path = Path(raw_job)
    if not job_path.is_absolute():
        parts = job_path.parts
        if not parts or parts[0] != "results":
            job_path = Path("results") / job_path
    job_dir = str(job_path.resolve())
    os.makedirs(job_dir, exist_ok=True)

    # Load base configuration if available
    cfg = {}
    if config_file and os.path.exists(config_file):
        try:
            with open(config_file, "r") as f:
                cfg = json.load(f)
        except Exception as e:
            print(f"[Notice] Could not load '{config_file}': {e}. Using defaults.")

    arm_cfg = cfg.get("arm_design", {})
    plp_cfg = cfg.get("plp_generation", {})

    target_count = int(count if count is not None else plp_cfg.get("plp_per_arm", 100))

    # Output file paths
    if not output_csv:
        final_csv_path = os.path.join(job_dir, "designed_plp_probes.csv")
    elif os.path.isabs(output_csv) or os.path.dirname(output_csv):
        final_csv_path = os.path.abspath(output_csv)
    else:
        final_csv_path = os.path.join(job_dir, output_csv)

    ranking_csv_path = os.path.join(job_dir, "shared_backbone_ranking.csv")
    combined_csv_path = os.path.join(job_dir, "combined_oligos.csv")

    configured_struct_dir = structures_dir or plp_cfg.get("structures_dir") or "plp_structures"
    if os.path.isabs(configured_struct_dir) or os.path.dirname(configured_struct_dir):
        final_struct_dir = os.path.abspath(configured_struct_dir)
    else:
        final_struct_dir = os.path.join(job_dir, configured_struct_dir)

    vienna_dir = os.path.join(final_struct_dir, "vienna")
    mfold_dir = os.path.join(final_struct_dir, "mfold")
    os.makedirs(vienna_dir, exist_ok=True)
    os.makedirs(mfold_dir, exist_ok=True)

    # Reference target (Target_1)
    target_1 = targets[0]
    t1_id = target_1["target_id"]
    clean_left = target_1["left_arm"]
    clean_right = target_1["right_arm"]
    avg_arm = (len(clean_left) + len(clean_right)) // 2

    # CRISPR resolution across targets
    pam_info_seed = {}
    clean_crrna_repeat = (crrna_repeat or "UAAUUUCUACUAAGUGUAGAU").strip().upper().replace("T", "U")
    if crispr_mode:
        pam_info_seed = determine_shared_pam_and_rest_tail(targets)

    # Architecture resolution
    num_threads = int(threads or plp_cfg.get("threads") or arm_cfg.get("threads") or os.cpu_count() or 32)
    arch_cfg = plp_cfg.get("backbone_architecture", {})
    s1_len = int(spacer1_len if spacer1_len is not None else arch_cfg.get("spacer1_len", 4))
    rca_len = int(rca_primer_binding_len if rca_primer_binding_len is not None else arch_cfg.get("rca_primer_binding_len", 20))
    s2_len = int(spacer2_len if spacer2_len is not None else arch_cfg.get("spacer2_len", 4))
    pr_len = int(probe_binding_len if probe_binding_len is not None else arch_cfg.get("probe_binding_len", 20))
    s3_len = int(spacer3_len if spacer3_len is not None else arch_cfg.get("spacer3_len", 0))
    r_len = int(rest_len if rest_len is not None else arch_cfg.get("rest_len", 4))
    univ_rca = (universal_rca_primer if universal_rca_primer is not None else arch_cfg.get("universal_rca_primer_seq", "")).strip()
    if univ_rca:
        rca_len = len(univ_rca)

    univ_probe = (universal_probe if universal_probe is not None else arch_cfg.get("universal_probe_seq", "")).strip()
    if univ_probe:
        pr_len = len(univ_probe)

    min_tm_filter = float(arm_cfg.get("min_tm", 40.0))
    re_site_str = str(restriction_site if restriction_site is not None else arch_cfg.get("restriction_site", "")).strip()
    re_seq = resolve_restriction_site(re_site_str) if re_site_str else ""
    re_len = len(re_seq)
    re_mode = arch_cfg.get("re_mode", "standalone")
    re_probe_offset = int(arch_cfg.get("re_probe_offset", 2))
    actual_re_len = re_len if re_mode == "standalone" else 0

    arch_bb_len = s1_len + rca_len + s2_len + pr_len + s3_len + actual_re_len + r_len

    eff_bb_factor = backbone_length_factor if backbone_length_factor is not None else plp_cfg.get("backbone_length_factor")
    if backbone_len is not None:
        target_bb_len = int(backbone_len)
    elif plp_cfg.get("backbone_length"):
        target_bb_len = int(plp_cfg.get("backbone_length"))
    elif eff_bb_factor is not None and str(eff_bb_factor).strip():
        target_bb_len = int(eff_bb_factor) * avg_arm
    else:
        target_bb_len = arch_bb_len

    if target_bb_len < arch_bb_len:
        actual_bb_len = arch_bb_len
    else:
        actual_bb_len = target_bb_len
        r_len += (actual_bb_len - arch_bb_len)
        arch_bb_len = actual_bb_len

    # Ensure r_len is at least the length of the engineered rest_tail in CRISPR mode
    if crispr_mode and pam_info_seed.get("rest_tail"):
        req_tail_len = len(pam_info_seed["rest_tail"])
        if r_len < req_tail_len:
            actual_bb_len += (req_tail_len - r_len)
            r_len = req_tail_len
            arch_bb_len = actual_bb_len

    parsed_screen_seqs = set()
    if screen_against:
        for inp in screen_against:
            if os.path.isfile(inp) and inp.lower().endswith(".csv"):
                with open(inp, "r") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        for col in ["Left_Arm", "Right_Arm", "PLP"]:
                            if col in row and row[col]:
                                parsed_screen_seqs.add(row[col].strip().upper())
            else:
                parsed_screen_seqs.add(inp.strip().upper())

    # Do not screen against reference arms
    parsed_screen_seqs.discard(clean_left)
    parsed_screen_seqs.discard(clean_right)

    generation_cfg = {
        "generate_plots": False,
        "screen_against": list(parsed_screen_seqs),
        "max_cross_dimer_dg": float(max_cross_dimer_dg) if max_cross_dimer_dg is not None else float(plp_cfg.get("max_cross_dimer_dg", -5.0)),
        "plp_per_arm": target_count,
        "threads": num_threads,
        "backbone_length": actual_bb_len,
        "maxHomopolymer": max(int(plp_cfg.get("maxHomopolymer", arm_cfg.get("maxHomopolymer", 4))), 4 if crispr_mode else 1),
        "max_hairpin_tm": float(plp_cfg.get("max_hairpin_tm", arm_cfg.get("max_hairpin_tm", 30.0))),
        "min_plp_mfe_dg": float(plp_cfg.get("min_plp_mfe_dg", -6.0)),
        "max_paired_arm_bases": int(plp_cfg.get("max_paired_arm_bases", 4)),
        "max_paired_total_bases": int(plp_cfg.get("max_paired_total_bases", 16)),
        "blastMinIdentity": float(plp_cfg.get("blastMinIdentity", arm_cfg.get("blastMinIdentity", 70.0))),
        "maxAllowedHits": plp_cfg.get("maxAllowedHits", arm_cfg.get("maxAllowedHits", 0)),
        "keep_ct_ps": keep_ct_ps or bool(plp_cfg.get("keep_ct_ps", False)),
        "backbone_architecture": {
            "spacer1_len": s1_len,
            "rca_primer_binding_len": rca_len,
            "spacer2_len": s2_len,
            "probe_binding_len": pr_len,
            "spacer3_len": s3_len,
            "rest_len": r_len,
            "restriction_site": re_site_str,
            "re_mode": re_mode,
            "re_probe_offset": re_probe_offset,
            "universal_rca_primer_seq": univ_rca,
            "universal_probe_seq": univ_probe,
            "min_rca_primer_tm": min_tm_filter,
            "min_probe_tm": min_tm_filter,
            "crispr_mode": crispr_mode,
            "rest_tail": pam_info_seed.get("rest_tail", "") if crispr_mode else "",
        },
    }

    blast_dbs = {}
    if not skip_blast and blastdb_dir and os.path.exists(blastdb_dir):
        blast_dbs = discover_blast_dbs(blastdb_dir)

    print("=" * 80)
    print("UNIVERSAL SHARED BACKBONE PADLOCK PROBE DESIGNER" + (" [CRISPR Cas12a Mode]" if crispr_mode else ""))
    print("=" * 80)
    print(f"Job Output Directory   : {job_dir}")
    print(f"Input Arms CSV         : {arms_csv} ({num_targets} targets detected)")
    print(f"Reference Seed Target  : {t1_id} (Left: {clean_left}, Right: {clean_right})")
    for idx, t in enumerate(targets[1:], start=2):
        print(f"Target {idx:<15}: {t['target_id']} (Left: {t['left_arm']}, Right: {t['right_arm']})")
    if crispr_mode:
        print(f"CRISPR Cas Mode        : ENABLED")
        print(f"Panel PAM Strategy     : {pam_info_seed.get('pam_seq', '[A/C/G]TT[A/C/G]')} ({pam_info_seed.get('pam_location', '')})")
        print(f"Engineered Rest Tail   : '{pam_info_seed.get('rest_tail', '')}' (appended at Rest 3'-end)")
        print(f"Cis-Probe Generation   : {'Enabled' if generate_cis_probe else 'Disabled (Default)'}")
        print(f"crRNA Repeat Sequence  : {clean_crrna_repeat} ({len(clean_crrna_repeat)} nt)")
        print(f"crRNA Spacer Length    : {spacer_len} nt")
        print(f"CRISPR FQ Reporter     : 5'-FAM-TTATTATT-BHQ1-3'")
        print(f"CRISPR FB Reporter     : 5'-FAM-TTATTATT-Biotin-3'")
    print(f"Candidate Backbones    : {target_count} to generate on reference {t1_id}")
    print(f"Backbone Architecture  : Spacer1({s1_len}nt) + RCA({rca_len}nt) + Spacer2({s2_len}nt) + Probe({pr_len}nt) + Spacer3({s3_len}nt) + RE({actual_re_len}nt) + Rest({r_len}nt) = {actual_bb_len} nt")
    print(f"Openness Criteria      : ViennaRNA dG >= {generation_cfg['min_plp_mfe_dg']} kcal/mol, Arm Paired <= {generation_cfg['max_paired_arm_bases']} nt, Total Paired <= {generation_cfg['max_paired_total_bases']} nt")
    print(f"Active CPU Cores       : {num_threads}")
    print(f"BLAST Background DBs   : {', '.join(blast_dbs.keys()) if blast_dbs else ('Skipped (--skip-blast)' if skip_blast else 'None found')}")
    print("-" * 80 + "\n")

    # Step 1: Generate candidate backbones for Target_1
    print(f"[Step 1/3] Generating {target_count} candidate open backbones for seed target '{t1_id}'...")
    total_counter = mp.Value("L", 0)
    lock = mp.Lock()
    stop_event = mp.Event()

    with mp.Pool(
        processes=num_threads,
        initializer=init_worker,
        initargs=(total_counter, lock, stop_event),
    ) as pool:
        records = generate_plp_for_arm(
            left_arm=clean_left,
            right_arm=clean_right,
            arm_meta={"Candidate_ID": f"{t1_id}_Cand"},
            cfg=generation_cfg,
            structures_dir=final_struct_dir,
            blast_dbs=blast_dbs,
            pool=pool,
            total_counter=total_counter,
            lock=lock,
            stop_event=stop_event,
        )

    print(f">> Successfully generated {len(records)} candidate backbones for {t1_id}.\n")

    # Step 2: Graft onto Target_2 ... Target_N and evaluate thermodynamics
    print(f"[Step 2/3] Grafting backbones across all {num_targets} targets and evaluating folding thermodynamics...")
    evaluated_backbones = []
    for cand_idx, cand in enumerate(records, start=1):
        bb = cand["bb"]
        target_evals = []
        for t_idx, t in enumerate(targets):
            t_id = t["target_id"]
            t_left = t["left_arm"]
            t_right = t["right_arm"]

            crispr_info = {}
            if crispr_mode:
                crispr_info = design_crispr_components(
                    left_arm=t_left,
                    right_arm=t_right,
                    backbone=bb,
                    crrna_repeat=clean_crrna_repeat,
                    spacer_len=spacer_len,
                    cis_probe_len=cis_probe_len,
                    generate_cis_probe=generate_cis_probe,
                )

            crispr_valid = crispr_info.get("valid_pam", True) if crispr_mode else True

            if t_idx == 0:
                db = cand["dot_bracket"]
                dg = cand["mfe_dg"]
                pl = cand["paired_left"]
                pr = cand["paired_right"]
                pt = cand["paired_total"]
                full_plp = cand["full_plp"]
                full_gc = cand["full_gc"]
                full_tm = cand["full_tm"]
                passes = cand.get("pass", True) and crispr_valid
            else:
                full_plp = t_left + bb + t_right
                db, dg = eval_plp_structure(full_plp, folding_engine="RNAfold")
                pl = sum(1 for ch in db[:len(t_left)] if ch in "()")
                pr = sum(1 for ch in db[-len(t_right):] if ch in "()")
                pt = sum(1 for ch in db if ch in "()")
                passes = (
                    dg >= generation_cfg["min_plp_mfe_dg"]
                    and pl <= generation_cfg["max_paired_arm_bases"]
                    and pr <= generation_cfg["max_paired_arm_bases"]
                    and pt <= generation_cfg["max_paired_total_bases"]
                    and crispr_valid
                )
                full_gc, full_tm = compute_metrics(Seq(full_plp))

            target_evals.append({
                "target_id": t_id,
                "left_arm": t_left,
                "right_arm": t_right,
                "full_plp": full_plp,
                "full_gc": full_gc,
                "full_tm": full_tm,
                "dot_bracket": db,
                "mfe_dg": dg,
                "paired_left": pl,
                "paired_right": pr,
                "paired_total": pt,
                "pass": passes,
                "vienna_plot": "",
                "crispr": crispr_info,
            })

        pass_cnt = sum(1 for e in target_evals if e["pass"])
        all_pass = (pass_cnt == num_targets)
        worst_dg = min(e["mfe_dg"] for e in target_evals)
        best_dg = max(e["mfe_dg"] for e in target_evals)
        avg_dg = sum(e["mfe_dg"] for e in target_evals) / float(num_targets)
        max_arm_p = max(max(e["paired_left"], e["paired_right"]) for e in target_evals)
        avg_arm_p = sum(e["paired_left"] + e["paired_right"] for e in target_evals) / (2.0 * float(num_targets))
        worst_total_p = max(e["paired_total"] for e in target_evals)

        evaluated_backbones.append({
            "cand_idx": cand_idx,
            "cand": cand,
            "bb": bb,
            "target_evals": target_evals,
            "pass_count": pass_cnt,
            "all_pass": all_pass,
            "worst_dg": worst_dg,
            "best_dg": best_dg,
            "avg_dg": avg_dg,
            "max_arm_paired": max_arm_p,
            "avg_arm_paired": avg_arm_p,
            "worst_total_paired": worst_total_p,
        })

    # Step 3: Composite Sorting and Ranking
    print(f"[Step 3/3] Ranking backbones across all {num_targets} targets...")
    # Preliminary sorting based on ViennaRNA openness metrics
    evaluated_backbones.sort(
        key=lambda x: (
            1 if x["all_pass"] else 0,
            x["pass_count"],
            x["worst_dg"],
            x["avg_dg"],
            -x["max_arm_paired"],
            -x["avg_arm_paired"],
            -x["worst_total_paired"],
        ),
        reverse=True,
    )

    all_pass_count = sum(1 for item in evaluated_backbones if item["all_pass"])
    print(f">> Found {all_pass_count}/{len(evaluated_backbones)} backbones that pass openness criteria on ALL {num_targets} targets!")

    # Evaluate DNA mfold for passing candidates (and top candidates) to perform fine-grained tie-breaking
    mfold_eval_candidates = [item for item in evaluated_backbones if item["all_pass"]]
    if len(mfold_eval_candidates) < 10:
        mfold_eval_candidates = evaluated_backbones[:min(10, len(evaluated_backbones))]

    print(f">> Evaluating DNA mfold thermodynamics for {len(mfold_eval_candidates)} candidate backbones for tie-breaking...")
    for item in mfold_eval_candidates:
        for te in item["target_evals"]:
            m_dg, m_db, _ = eval_structure_mfold_detailed(
                te["full_plp"],
                seq_id=f"tie_{item['cand_idx']}_{te['target_id']}",
                output_dir=None,
                temperature=37.0,
                keep_ct_ps=False,
            )
            m_pl = sum(1 for ch in m_db[:len(te["left_arm"])] if ch in "()")
            m_pr = sum(1 for ch in m_db[-len(te["right_arm"]):] if ch in "()")
            m_pt = sum(1 for ch in m_db if ch in "()")
            te["mfold_dg"] = m_dg if m_dg is not None else 0.0
            te["mfold_db"] = m_db
            te["mfold_paired_left"] = m_pl
            te["mfold_paired_right"] = m_pr
            te["mfold_paired_total"] = m_pt

        item["mfold_worst_dg"] = min(e["mfold_dg"] for e in item["target_evals"])
        item["mfold_avg_dg"] = sum(e["mfold_dg"] for e in item["target_evals"]) / float(num_targets)
        item["mfold_max_arm_paired"] = max(max(e["mfold_paired_left"], e["mfold_paired_right"]) for e in item["target_evals"])
        item["mfold_worst_total_paired"] = max(e["mfold_paired_total"] for e in item["target_evals"])

    # Default values for candidates not evaluated by mfold
    for item in evaluated_backbones:
        if "mfold_worst_dg" not in item:
            item["mfold_worst_dg"] = -999.0
            item["mfold_avg_dg"] = -999.0
            item["mfold_max_arm_paired"] = 999
            item["mfold_worst_total_paired"] = 999
            for te in item["target_evals"]:
                te["mfold_dg"] = 0.0
                te["mfold_db"] = ""
                te["mfold_paired_left"] = 0
                te["mfold_paired_right"] = 0
                te["mfold_paired_total"] = 0

    # Hierarchical re-sorting: Vienna criteria first, then DNA mfold tie-breaking:
    # 1. all_pass (True > False)
    # 2. pass_count (e.g. 2/2 > 1/2)
    # 3. Vienna worst_dg (closer to 0.0 kcal/mol is more open)
    # 4. Vienna avg_dg
    # 5. Vienna max_arm_paired (lower is better)
    # 6. Vienna avg_arm_paired (lower is better)
    # 7. Vienna worst_total_paired (lower is better)
    # 8. mfold max_arm_paired (lower is better, e.g. 0 bp > 1 bp)
    # 9. mfold worst_total_paired (lower is better)
    # 10. mfold worst_dg (closer to 0.0 or positive is more open)
    # 11. mfold avg_dg
    evaluated_backbones.sort(
        key=lambda x: (
            1 if x["all_pass"] else 0,
            x["pass_count"],
            x["worst_dg"],
            x["avg_dg"],
            -x["max_arm_paired"],
            -x["avg_arm_paired"],
            -x["worst_total_paired"],
            -x["mfold_max_arm_paired"],
            -x["mfold_worst_total_paired"],
            x["mfold_worst_dg"],
            x["mfold_avg_dg"],
        ),
        reverse=True,
    )

    # Clean up vienna & mfold directories so only Top 10 plots are preserved
    if os.path.exists(vienna_dir):
        for f in os.listdir(vienna_dir):
            if f.endswith(".svg"):
                try:
                    os.remove(os.path.join(vienna_dir, f))
                except OSError:
                    pass
    if os.path.exists(mfold_dir):
        for f in os.listdir(mfold_dir):
            try:
                os.remove(os.path.join(mfold_dir, f))
            except OSError:
                pass

    # Assign ranks and render Vienna SVG plots & DNA mfold PNG plots strictly for Top 10 backbones across all targets
    top_plot_n = min(10, len(evaluated_backbones))
    for rank, item in enumerate(evaluated_backbones, start=1):
        item["rank"] = rank
        bb_id = f"Shared_BB_{rank}"
        item["bb_id"] = bb_id

        if rank <= top_plot_n:
            for te in item["target_evals"]:
                plot_id = f"{te['target_id']}_Rank_{rank}"
                svg_path = render_vienna_svg(plot_id, te["full_plp"], te["dot_bracket"], vienna_dir)
                te["vienna_plot"] = svg_path

                m_dg, m_db, m_plot = eval_structure_mfold_detailed(
                    te["full_plp"],
                    seq_id=plot_id,
                    output_dir=mfold_dir,
                    temperature=37.0,
                    keep_ct_ps=keep_ct_ps,
                )
                te["mfold_dg"] = m_dg if m_dg is not None else 0.0
                te["mfold_db"] = m_db
                te["mfold_plot"] = m_plot
                te["mfold_paired_left"] = sum(1 for ch in m_db[:len(te["left_arm"])] if ch in "()")
                te["mfold_paired_right"] = sum(1 for ch in m_db[-len(te["right_arm"]):] if ch in "()")
                te["mfold_paired_total"] = sum(1 for ch in m_db if ch in "()")
        else:
            for te in item["target_evals"]:
                te["vienna_plot"] = ""
                te["mfold_plot"] = ""

    print(f">> Generated Vienna SVG and DNA mfold structure plots for Top {top_plot_n} backbones across all {num_targets} targets.\n")

    # Export CSV 1: shared_backbone_ranking.csv
    ranking_headers = [
        "Rank",
        "Backbone_ID",
        "All_Targets_Pass",
        "Pass_Count",
        "Worst_dG_kcal_mol",
        "Avg_dG_kcal_mol",
        "Best_dG_kcal_mol",
        "Max_Arm_Paired_Bases",
        "Worst_Total_Paired_Bases",
        "mfold_Worst_dG_kcal_mol",
        "mfold_Max_Arm_Paired_Bases",
        "mfold_Worst_Total_Paired_Bases",
    ]
    for t in targets:
        tid = t["target_id"]
        ranking_headers.extend([
            f"{tid}_dG_kcal_mol",
            f"{tid}_Pass",
            f"{tid}_Paired_Left",
            f"{tid}_Paired_Right",
            f"{tid}_Paired_Total",
            f"{tid}_Dot_Bracket",
            f"{tid}_mfold_dG_kcal_mol",
            f"{tid}_mfold_Paired_Left",
            f"{tid}_mfold_Paired_Right",
            f"{tid}_mfold_Paired_Total",
            f"{tid}_mfold_Dot_Bracket",
        ])
        if crispr_mode:
            ranking_headers.append(f"{tid}_PAM")
            if generate_cis_probe:
                ranking_headers.append(f"{tid}_Cis_Probe")
            ranking_headers.append(f"{tid}_crRNA_Spacer")
    ranking_headers.extend([
        "Backbone_Sequence",
        "Backbone_Length",
        "Backbone_GC%",
        "Backbone_Tm",
        "Backbone_Max_Homopolymer",
        "Spacer1_Seq",
        "RCA_Primer_Binding_Region",
        "Spacer2_Seq",
        "Probe_Binding_Region",
        "Rest_Backbone_Seq",
        "RCA_Primer_Seq_5to3",
        "RCA_Primer_5_Mod",
        "RCA_Primer_Tm",
        "RCA_Primer_GC%",
        "Detection_Probe_Seq_5to3",
        "Detection_Probe_5_Mod",
        "Detection_Probe_Tm",
        "Detection_Probe_GC%",
    ])
    for db_name in blast_dbs:
        ranking_headers.append(f"Backbone_Hits_{db_name}")

    with open(ranking_csv_path, "w", newline="") as f_rank:
        writer = csv.writer(f_rank)
        writer.writerow(ranking_headers)
        for item in evaluated_backbones:
            cand = item["cand"]
            m_wdg = item.get("mfold_worst_dg", -999.0)
            m_map = item.get("mfold_max_arm_paired", 999)
            m_wtp = item.get("mfold_worst_total_paired", 999)
            row = [
                item["rank"],
                item["bb_id"],
                item["all_pass"],
                f"{item['pass_count']}/{num_targets}",
                f"{item['worst_dg']:.2f}",
                f"{item['avg_dg']:.2f}",
                f"{item['best_dg']:.2f}",
                item["max_arm_paired"],
                item["worst_total_paired"],
                f"{m_wdg:.2f}" if m_wdg != -999.0 else "N/A",
                m_map if m_map != 999 else "N/A",
                m_wtp if m_wtp != 999 else "N/A",
            ]
            for te in item["target_evals"]:
                row.extend([
                    f"{te['mfe_dg']:.2f}",
                    te["pass"],
                    te["paired_left"],
                    te["paired_right"],
                    te["paired_total"],
                    te["dot_bracket"],
                    f"{te.get('mfold_dg', 0.0):.2f}" if te.get("mfold_db") else "N/A",
                    te.get("mfold_paired_left", "N/A") if te.get("mfold_db") else "N/A",
                    te.get("mfold_paired_right", "N/A") if te.get("mfold_db") else "N/A",
                    te.get("mfold_paired_total", "N/A") if te.get("mfold_db") else "N/A",
                    te.get("mfold_db", "N/A") if te.get("mfold_db") else "N/A",
                ])
                if crispr_mode:
                    c = te.get("crispr", {})
                    row.append(c.get("pam_seq", ""))
                    if generate_cis_probe:
                        row.append(c.get("cis_probe_seq", ""))
                    row.append(c.get("spacer_seq_rna", ""))
            row.extend([
                item["bb"],
                len(item["bb"]),
                cand.get("bb_gc", ""),
                cand.get("bb_tm", ""),
                cand.get("bb_homo", ""),
                cand.get("Spacer1_Seq", ""),
                cand.get("RCA_Primer_Binding_Region", ""),
                cand.get("Spacer2_Seq", ""),
                cand.get("Probe_Binding_Region", ""),
                cand.get("Rest_Backbone_Seq", ""),
                cand.get("RCA_Primer_Seq_5to3", ""),
                cand.get("RCA_Primer_5_Mod", "Biotin"),
                cand.get("RCA_Primer_Tm", ""),
                cand.get("RCA_Primer_GC%", ""),
                cand.get("Detection_Probe_Seq_5to3", ""),
                cand.get("Detection_Probe_5_Mod", "6-FAM"),
                cand.get("Detection_Probe_Tm", ""),
                cand.get("Detection_Probe_GC%", ""),
            ])
            row.extend(cand.get("blast_hits", []))
            writer.writerow(row)

    # Export CSV 2: designed_plp_probes.csv
    probe_headers = [
        "Rank",
        "Backbone_ID",
        "Target_ID",
        "PLP_ID",
        "Full_PLP_Sequence",
        "Total_Length",
        "Full_PLP_GC%",
        "Full_PLP_Tm",
        "Left_Arm",
        "Backbone",
        "Right_Arm",
        "Backbone_Length",
        "Backbone_GC%",
        "Backbone_Tm",
        "Spacer1_Seq",
        "RCA_Primer_Binding_Region",
        "Spacer2_Seq",
        "Probe_Binding_Region",
        "Rest_Backbone_Seq",
        "RCA_Primer_Seq_5to3",
        "RCA_Primer_5_Mod",
        "RCA_Primer_Tm",
        "RCA_Primer_GC%",
        "Detection_Probe_Seq_5to3",
        "Detection_Probe_5_Mod",
        "Detection_Probe_Tm",
        "Detection_Probe_GC%",
        "Vienna_Full_PLP_Structure",
        "Vienna_Full_PLP_dG_kcal_mol",
        "mfold_Full_PLP_Structure",
        "mfold_Full_PLP_dG_kcal_mol",
        "Left_Arm_Paired_Bases",
        "Right_Arm_Paired_Bases",
        "Total_Paired_Bases",
        "mfold_Left_Arm_Paired_Bases",
        "mfold_Right_Arm_Paired_Bases",
        "mfold_Total_Paired_Bases",
        "Pass_Openness",
        "All_Targets_Pass",
        "Vienna_Structure_Plot",
        "mfold_Structure_Plot",
    ]
    if crispr_mode:
        probe_headers.extend([
            "CRISPR_PAM_Motif",
            "CRISPR_PAM_Location",
        ])
        if generate_cis_probe:
            probe_headers.extend([
                "CRISPR_Cis_Probe_Seq_5to3",
                "CRISPR_Cis_Probe_Length",
                "CRISPR_Cis_Probe_Tm",
                "CRISPR_Cis_Probe_GC%",
            ])
        probe_headers.extend([
            "CRISPR_crRNA_Full_Seq_5to3",
            "CRISPR_crRNA_Repeat",
            "CRISPR_crRNA_Spacer",
            "CRISPR_crRNA_Spacer_Length",
        ])
    with open(final_csv_path, "w", newline="") as f_probes:
        writer = csv.writer(f_probes)
        writer.writerow(probe_headers)
        for item in evaluated_backbones:
            cand = item["cand"]
            for te in item["target_evals"]:
                plp_id = f"{te['target_id']}_Rank{item['rank']}"
                row = [
                    item["rank"],
                    item["bb_id"],
                    te["target_id"],
                    plp_id,
                    te["full_plp"],
                    len(te["full_plp"]),
                    te["full_gc"],
                    te["full_tm"],
                    te["left_arm"],
                    item["bb"],
                    te["right_arm"],
                    len(item["bb"]),
                    cand.get("bb_gc", ""),
                    cand.get("bb_tm", ""),
                    cand.get("Spacer1_Seq", ""),
                    cand.get("RCA_Primer_Binding_Region", ""),
                    cand.get("Spacer2_Seq", ""),
                    cand.get("Probe_Binding_Region", ""),
                    cand.get("Rest_Backbone_Seq", ""),
                    cand.get("RCA_Primer_Seq_5to3", ""),
                    cand.get("RCA_Primer_5_Mod", "Biotin"),
                    cand.get("RCA_Primer_Tm", ""),
                    cand.get("RCA_Primer_GC%", ""),
                    cand.get("Detection_Probe_Seq_5to3", ""),
                    cand.get("Detection_Probe_5_Mod", "6-FAM"),
                    cand.get("Detection_Probe_Tm", ""),
                    cand.get("Detection_Probe_GC%", ""),
                    te["dot_bracket"],
                    f"{te['mfe_dg']:.2f}",
                    te.get("mfold_db", ""),
                    f"{te.get('mfold_dg', 0.0):.2f}",
                    te["paired_left"],
                    te["paired_right"],
                    te["paired_total"],
                    te.get("mfold_paired_left", "N/A"),
                    te.get("mfold_paired_right", "N/A"),
                    te.get("mfold_paired_total", "N/A"),
                    te["pass"],
                    item["all_pass"],
                    te.get("vienna_plot", ""),
                    te.get("mfold_plot", ""),
                ]
                if crispr_mode:
                    c = te.get("crispr", {})
                    row.extend([
                        c.get("pam_seq", ""),
                        c.get("pam_location", ""),
                    ])
                    if generate_cis_probe:
                        row.extend([
                            c.get("cis_probe_seq", ""),
                            c.get("cis_probe_len", ""),
                            c.get("cis_probe_tm", ""),
                            c.get("cis_probe_gc", ""),
                        ])
                    row.extend([
                        c.get("crrna_full_seq", ""),
                        c.get("crrna_repeat", ""),
                        c.get("spacer_seq_rna", ""),
                        c.get("spacer_len", ""),
                    ])
                writer.writerow(row)

    # Export CSV 3: combined_oligos.csv (Unified All-in-One Assay Oligo Sheet)
    combined_csv_path = os.path.join(job_dir, "combined_oligos.csv")
    combined_headers = [
        "Oligo_Name",
        "Sequence_5to3",
        "Modifications",
        "Type",
        "Length_nt",
        "GC%",
        "Tm_C",
        "Assay_Category",
        "Role",
        "Associated_Rank",
        "Associated_Backbone_ID",
        "Target_ID",
    ]
    with open(combined_csv_path, "w", newline="") as f_comb:
        writer = csv.writer(f_comb)
        writer.writerow(combined_headers)
        for item in evaluated_backbones[:3]:
            rank = item["rank"]
            bb_id = item["bb_id"]
            cand = item["cand"]
            rank_label = f"Rank {rank}" + (" (Primary Selection)" if rank == 1 else " (Backup)")

            # 1. Target PLPs
            for te in item["target_evals"]:
                tid = te["target_id"]
                writer.writerow([
                    f"{tid}_PLP_{bb_id}",
                    te["full_plp"],
                    "5'-Phos",
                    "DNA",
                    len(te["full_plp"]),
                    te["full_gc"],
                    te["full_tm"],
                    "PLP Ligation Probe",
                    f"Padlock Probe for {tid} [{rank_label}]",
                    rank,
                    bb_id,
                    tid,
                ])

            # 2. Universal RCA Primer
            if cand.get("RCA_Primer_Seq_5to3"):
                writer.writerow([
                    f"Universal_RCA_Primer_{bb_id}",
                    cand["RCA_Primer_Seq_5to3"],
                    "5'-Biotin",
                    "DNA",
                    len(cand["RCA_Primer_Seq_5to3"]),
                    cand.get("RCA_Primer_GC%", ""),
                    cand.get("RCA_Primer_Tm", ""),
                    "Universal RCA Primer",
                    f"Universal RCA Amplification Primer (Shared) [{rank_label}]",
                    rank,
                    bb_id,
                    "Universal (All Targets)",
                ])

            # 3. Universal Detection Probe
            if cand.get("Detection_Probe_Seq_5to3"):
                writer.writerow([
                    f"Universal_Detection_Probe_{bb_id}",
                    cand["Detection_Probe_Seq_5to3"],
                    "5'-6-FAM",
                    "DNA",
                    len(cand["Detection_Probe_Seq_5to3"]),
                    cand.get("Detection_Probe_GC%", ""),
                    cand.get("Detection_Probe_Tm", ""),
                    "Universal LFA FAM Probe",
                    f"Universal LFA Fluorescent Probe (Shared) [{rank_label}]",
                    rank,
                    bb_id,
                    "Universal (All Targets)",
                ])

            # 4. If crispr_mode: Target crRNA Guides & Cis Probes
            if crispr_mode:
                for te in item["target_evals"]:
                    tid = te["target_id"]
                    c = te.get("crispr", {})
                    if generate_cis_probe and c.get("cis_probe_seq"):
                        writer.writerow([
                            f"{tid}_Cis_Probe_{bb_id}",
                            c.get("cis_probe_seq", ""),
                            "None",
                            "DNA",
                            c.get("cis_probe_len", len(c.get("cis_probe_seq", ""))),
                            c.get("cis_probe_gc", ""),
                            c.get("cis_probe_tm", ""),
                            "CRISPR Cis Detection Probe",
                            f"CRISPR Cas12a Activator Cis-Probe for {tid} ({c.get('pam_location', '')}) [{rank_label}]",
                            rank,
                            bb_id,
                            tid,
                        ])
                    if c.get("crrna_full_seq"):
                        writer.writerow([
                            f"{tid}_crRNA_{bb_id}",
                            c.get("crrna_full_seq"),
                            "None (RNA)",
                            "RNA",
                            len(c.get("crrna_full_seq", "")),
                            c.get("spacer_gc", ""),
                            c.get("spacer_tm", ""),
                            "CRISPR Cas12a Guide RNA",
                            f"crRNA Guide for {tid} (Repeat: {len(c.get('crrna_repeat', ''))}nt, Spacer: {c.get('spacer_len', '')}nt) [{rank_label}]",
                            rank,
                            bb_id,
                            tid,
                        ])

        # 5. Universal CRISPR Reporters
        if crispr_mode:
            for rep in CRISPR_REPORTERS:
                cat = "CRISPR Fluorescent Reporter (FQ)" if "FQ" in rep["oligo_name"] else "CRISPR Lateral Flow Reporter (FB)"
                writer.writerow([
                    rep["oligo_name"],
                    rep["sequence_5to3"],
                    rep["modifications"],
                    rep["type"],
                    rep["length_nt"],
                    rep["gc_percent"],
                    rep["tm_c"],
                    cat,
                    rep["role"],
                    "Universal",
                    "Shared",
                    "Universal (All Targets)",
                ])

    # Display Top 10 Summary Table
    top_n = min(10, len(evaluated_backbones))
    print(f"Top {top_n} Universal Shared Backbones (Ranked across all {num_targets} targets):")
    hdr = f"{'Rank':<6} {'Pass Rate':<11} {'Vienna dG':<12} {'mfold dG':<12} {'Vienna Arm':<12} {'mfold Arm':<11}"
    for t in targets:
        hdr += f" {t['target_id'][:15] + ' dG':<18}"
    hdr += f" {'RCA Primer (Tm)':<20} {'FAM Probe (Tm)':<20}"
    print(hdr)
    print("-" * len(hdr))

    for item in evaluated_backbones[:top_n]:
        cand = item["cand"]
        pass_str = f"{item['pass_count']}/{num_targets} Pass"
        m_wdg_val = item.get("mfold_worst_dg", -999.0)
        m_wdg_str = f"{m_wdg_val:>6.2f} kcal" if m_wdg_val != -999.0 else "    N/A    "
        m_arm_val = item.get("mfold_max_arm_paired", 999)
        m_arm_str = f"{m_arm_val:>3} bp" if m_arm_val != 999 else " N/A "
        row_str = f"#{item['rank']:<5} {pass_str:<11} {item['worst_dg']:>6.2f} kcal  {m_wdg_str:<12} {item['max_arm_paired']:>5} bp     {m_arm_str:<11}"
        for te in item["target_evals"]:
            row_str += f" {te['mfe_dg']:>8.2f} kcal     "
        rca_str = f"{cand.get('RCA_Primer_Seq_5to3', '')[:6]}.. ({cand.get('RCA_Primer_Tm', '')}°C)"
        probe_str = f"{cand.get('Detection_Probe_Seq_5to3', '')[:6]}.. ({cand.get('Detection_Probe_Tm', '')}°C)"
        row_str += f" {rca_str:<20} {probe_str:<20}"
        print(row_str)
    print("-" * len(hdr) + "\n")

    # Detailed report of Rank #1 Universal Shared Backbone
    best = evaluated_backbones[0]
    best_cand = best["cand"]
    print("=" * 80)
    print("RANK #1 WINNING UNIVERSAL SHARED BACKBONE DETAILS")
    print("=" * 80)
    print(f"Backbone ID              : {best['bb_id']}")
    print(f"Backbone Sequence (5'->3'): {best['bb']}")
    print(f"Backbone Length          : {len(best['bb'])} nt (GC: {best_cand.get('bb_gc')}%, Tm: {best_cand.get('bb_tm')}°C)")
    print(f"Pass Rate Across Targets : {best['pass_count']}/{num_targets} Targets Pass Openness Criteria")
    print(f"Thermodynamic Profile    : Worst dG: {best['worst_dg']:.2f} kcal/mol | Avg dG: {best['avg_dg']:.2f} kcal/mol")
    print(f"Max Paired Arm Bases     : {best['max_arm_paired']} bp (Arms remain accessible for hybridizing)")
    print(f"Universal RCA Primer     : {best_cand.get('RCA_Primer_Seq_5to3')} ({len(best_cand.get('RCA_Primer_Seq_5to3', ''))} nt, Tm: {best_cand.get('RCA_Primer_Tm')}°C, 5'-Biotin)")
    print(f"Universal FAM Probe      : {best_cand.get('Detection_Probe_Seq_5to3')} ({len(best_cand.get('Detection_Probe_Seq_5to3', ''))} nt, Tm: {best_cand.get('Detection_Probe_Tm')}°C, 5'-6-FAM)")
    print("\nTarget-Specific Padlock Probes (Ready to synthesize):")
    for te in best["target_evals"]:
        print(f"  * {te['target_id']} PLP ({len(te['full_plp'])} nt, Vienna dG: {te['mfe_dg']:.2f} kcal/mol, Paired Left: {te['paired_left']} bp, Paired Right: {te['paired_right']} bp, Status: {'PASS' if te['pass'] else 'FAIL'}):")
        print(f"    5'-{te['full_plp']}-3'")

    if crispr_mode:
        print("\nCRISPR CAS12A ASSAY OLIGOS (RANK #1 SELECTION ACROSS ALL TARGETS):")
        if generate_cis_probe:
            print(f"{'Target ID':<18} {'PAM':<6} {'Cis Probe (5\'->3\')':<34} {'crRNA Spacer (5\'->3\')':<24}")
            print("-" * 86)
            for te in best["target_evals"]:
                c = te.get("crispr", {})
                print(f"{te['target_id']:<18} {c.get('pam_seq', ''):<6} {c.get('cis_probe_seq', ''):<34} {c.get('spacer_seq_rna', ''):<24}")
            print("-" * 86)
        else:
            print(f"{'Target ID':<18} {'PAM':<6} {'crRNA Spacer (5\'->3\')':<24} {'PAM Location':<30}")
            print("-" * 82)
            for te in best["target_evals"]:
                c = te.get("crispr", {})
                print(f"{te['target_id']:<18} {c.get('pam_seq', ''):<6} {c.get('spacer_seq_rna', ''):<24} {c.get('pam_location', ''):<30}")
            print("-" * 82)
        print("Universal CRISPR Cas12a Reporter Probes:")
        for rep in CRISPR_REPORTERS:
            print(f"  * {rep['oligo_name']:<26}: {rep['sequence_5to3']} ({rep['role']})")

    print("=" * 80)

    print("\n" + "=" * 80)
    print("SUCCESS: Designed & Ranked Shared Universal Backbone across all targets!")
    print(f"Job Output Directory   -> {job_dir}")
    print(f"Backbone Ranking CSV   -> {ranking_csv_path}")
    print(f"Designed Probes CSV    -> {final_csv_path}")
    print(f"Combined Oligos Sheet  -> {combined_csv_path}")
    print(f"Vienna SVG Plots       -> {os.path.abspath(vienna_dir)}/")
    print(f"DNA mfold PNG Plots    -> {os.path.abspath(mfold_dir)}/")
    print("=" * 80 + "\n")

    return evaluated_backbones


def main():
    parser = argparse.ArgumentParser(
        description="Standalone Padlock Probe (PLP) Designer: Generate open, non-cross-reactive PLPs from custom Left and Right Arms (or multiple targets via CSV) into a designated job directory."
    )
    parser.add_argument("-j", "--job-name", required=True, help="Job name / output directory where probes and structures will be saved")
    parser.add_argument("-l", "--left-arm", default=None, help="Left hybridizing arm sequence (5' -> 3'). Required if --arms-csv is not provided.")
    parser.add_argument("-r", "--right-arm", default=None, help="Right hybridizing arm sequence (5' -> 3'). Required if --arms-csv is not provided.")
    parser.add_argument("-a", "--arms-csv", "--targets-csv", default=None, help="Input CSV containing multiple targets (Left_Arm and Right_Arm) to design a universal/shared backbone for.")
    parser.add_argument("-n", "--count", type=int, default=None, help="Number of open PLP probes/candidate backbones to generate (default: 100 for --arms-csv, or plp_per_arm/5 for single arm)")
    parser.add_argument(
        "-b",
        "--backbone-len",
        type=int,
        default=None,
        help="Custom backbone length in nt (default: computed from backbone_length_factor)",
    )
    parser.add_argument(
        "--backbone-length-factor",
        type=int,
        default=None,
        help="Multiplier for backbone size relative to average arm length (default: from config plp_generation.backbone_length_factor or 3)",
    )
    parser.add_argument(
        "-o",
        "--output-csv",
        default=None,
        help="Custom CSV filename inside job directory (default: designed_plp_probes.csv)",
    )
    parser.add_argument(
        "-s",
        "--structures-dir",
        default=None,
        help="Custom directory name for plots inside job directory (default: plp_structures)",
    )
    parser.add_argument(
        "--blastdb",
        default="blastdb",
        help="Directory containing background BLAST databases (default: blastdb)",
    )
    parser.add_argument(
        "--skip-blast",
        action="store_true",
        help="Skip BLAST background cross-reactivity screening (thermodynamic design only)",
    )
    parser.add_argument(
        "-c",
        "--config",
        default="config.json",
        help="Optional config.json to inherit thresholds from (default: config.json)",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=None,
        help="Number of CPU worker cores to allocate (default: all available cores)",
    )
    parser.add_argument(
        "--keep-ct-ps",
        action="store_true",
        help="Preserve raw .ct and .ps structure files from mfold",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        default="Custom_PLP",
        help="Probe ID prefix (default: Custom_PLP)",
    )
    parser.add_argument(
        "--universal-rca-primer",
        default=None,
        help="Optional fixed universal 5'->3' RCA primer sequence (5'-Biotin labeled)",
    )
    parser.add_argument(
        "--universal-probe",
        default=None,
        help="Optional fixed universal 5'->3' Detection probe sequence (5'-FAM labeled)",
    )
    parser.add_argument(
        "--screen-against",
        nargs="*",
        default=[],
        help="Sequences or CSV files containing arms/PLPs to screen against to prevent cross-reactivity",
    )
    parser.add_argument(
        "--max-cross-dimer-dg",
        type=float,
        default=None,
        help="Max allowed cross-reactivity dG (kcal/mol) with screen-against sequences",
    )
    parser.add_argument(
        "--spacer1-len",
        type=int,
        default=None,
        help="Length of Spacer 1 between Left Arm and RCA Primer Binding Region (default: 4 nt)",
    )
    parser.add_argument(
        "--spacer2-len",
        type=int,
        default=None,
        help="Length of Spacer 2 between RCA Primer and Probe Binding Region (default: 4 nt)",
    )
    parser.add_argument(
        "--rca-binding-len",
        type=int,
        default=None,
        help="Length of RCA Primer Binding Region (default: 20 nt)",
    )
    parser.add_argument(
        "--probe-binding-len",
        type=int,
        default=None,
        help="Length of Detection Probe Binding Region (default: 20 nt)",
    )
    parser.add_argument(
        "--rest-len",
        type=int,
        default=None,
        help="Length of Rest sequence between Probe Binding Region and Right Arm (default: 4 nt)",
    )
    parser.add_argument(
        "--restriction-site",
        default=None,
        help="Optional restriction enzyme cleavage site (e.g. BstUI, EcoRV, HaeIII, or raw sequence like CGCG) inserted into PLP backbone",
    )
    parser.add_argument(
        "--crispr",
        action="store_true",
        help="Enable CRISPR mode: engineers [A/C/G]-T-T-[A/C/G] (VTTV) PAM at Rest-Right Arm junction, designs crRNA guides, and exports reporter probes (5'-FAM-TTATTATT-BHQ1-3' and 5'-FAM-TTATTATT-Biotin-3').",
    )
    parser.add_argument(
        "--cis-probe",
        action="store_true",
        default=False,
        help="Enable CRISPR activator Cis-probe DNA oligo generation (default: False, cis-probe generation is off by default)",
    )
    parser.add_argument(
        "--crrna-repeat",
        default="UAAUUUCUACUAAGUGUAGAU",
        help="crRNA scaffold repeat RNA sequence (default: UAAUUUCUACUAAGUGUAGAU)",
    )
    parser.add_argument(
        "--spacer-len",
        type=int,
        default=20,
        help="crRNA guide spacer length in nt (default: 20 nt)",
    )
    parser.add_argument(
        "--cis-probe-len",
        type=int,
        default=28,
        help="CRISPR activator Cis probe target length in nt (default: 28 nt, 25-31 nt span; used when --cis-probe is enabled)",
    )

    args = parser.parse_args()

    crispr_active = bool(args.crispr or args.cis_probe)
    generate_cis = bool(args.cis_probe)

    try:
        if args.arms_csv:
            count = args.count if args.count is not None else 100
            design_shared_backbone_from_csv(
                arms_csv=args.arms_csv,
                job_name=args.job_name,
                count=count,
                backbone_len=args.backbone_len,
                backbone_length_factor=args.backbone_length_factor,
                output_csv=args.output_csv,
                structures_dir=args.structures_dir,
                blastdb_dir=args.blastdb,
                skip_blast=args.skip_blast,
                threads=args.threads,
                config_file=args.config,
                keep_ct_ps=args.keep_ct_ps,
                probe_prefix=args.prefix,
                universal_rca_primer=args.universal_rca_primer,
                universal_probe=args.universal_probe,
                spacer1_len=args.spacer1_len,
                spacer2_len=args.spacer2_len,
                rca_primer_binding_len=args.rca_binding_len,
                probe_binding_len=args.probe_binding_len,
                rest_len=args.rest_len,
                restriction_site=args.restriction_site,
                screen_against=args.screen_against,
                max_cross_dimer_dg=args.max_cross_dimer_dg,
                crispr_mode=crispr_active,
                crrna_repeat=args.crrna_repeat,
                spacer_len=args.spacer_len,
                cis_probe_len=args.cis_probe_len,
                generate_cis_probe=generate_cis,
            )
        else:
            if not args.left_arm or not args.right_arm:
                parser.error("Either --arms-csv (-a) or both --left-arm (-l) and --right-arm (-r) must be specified.")
            design_custom_plp(
                left_arm=args.left_arm,
                right_arm=args.right_arm,
                job_name=args.job_name,
                count=args.count,
                backbone_len=args.backbone_len,
                backbone_length_factor=args.backbone_length_factor,
                output_csv=args.output_csv,
                structures_dir=args.structures_dir,
                blastdb_dir=args.blastdb,
                skip_blast=args.skip_blast,
                threads=args.threads,
                config_file=args.config,
                keep_ct_ps=args.keep_ct_ps,
                probe_prefix=args.prefix,
                universal_rca_primer=args.universal_rca_primer,
                universal_probe=args.universal_probe,
                spacer1_len=args.spacer1_len,
                spacer2_len=args.spacer2_len,
                rca_primer_binding_len=args.rca_binding_len,
                probe_binding_len=args.probe_binding_len,
                rest_len=args.rest_len,
                restriction_site=args.restriction_site,
                screen_against=args.screen_against,
                max_cross_dimer_dg=args.max_cross_dimer_dg,
                crispr_mode=crispr_active,
                crrna_repeat=args.crrna_repeat,
                spacer_len=args.spacer_len,
                cis_probe_len=args.cis_probe_len,
                generate_cis_probe=generate_cis,
            )
    except Exception as err:
        print(f"\n[Error] {err}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
