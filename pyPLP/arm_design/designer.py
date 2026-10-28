import csv
import json
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq

from pyPLP.arm_design.aln_parser import parse_aln_pos_offset, parse_conservation_metadata
from pyPLP.arm_design.dimer_eval import eval_arm_dimer_vienna, eval_arm_dimer_mfold
from pyPLP.arm_design.junction_filter import evaluate_ligation_junction
from pyPLP.blast.discovery import discover_blast_dbs
from pyPLP.blast.search import run_blastn_short_parallel
from pyPLP.thermodynamics.hairpin import check_hairpin
from pyPLP.thermodynamics.homopolymer import get_max_homopolymer
from pyPLP.thermodynamics.metrics import compute_metrics


def design_plp_arms(input_fasta: str, output_csv: str, cfg: dict) -> bool:
    """Designs PLP candidate arms from conserved FASTA sequences.

    Applies thermodynamic (ViennaRNA), structural, and BLAST specificity filters,
    and profiles mfold dG on surviving candidates post-BLAST.
    """
    arm_cfg = cfg.get("arm_design", cfg)

    raw_lengths = arm_cfg.get("armLengths", [15])
    if isinstance(raw_lengths, int):
        arm_lengths = [raw_lengths]
    elif isinstance(raw_lengths, (list, tuple)):
        arm_lengths = [int(x) for x in raw_lengths]
    else:
        arm_lengths = [15]
    filter_by_tm = bool(arm_cfg.get("filter_by_tm", True))
    filter_by_homopolymer = bool(arm_cfg.get("filter_by_homopolymer", True))
    filter_by_hairpin = bool(arm_cfg.get("filter_by_hairpin", True))
    filter_by_arm_dimer = bool(arm_cfg.get("filter_by_arm_dimer", True))

    # Option A: Ligation Junction Quality Filter
    junc_cfg = arm_cfg.get("ligation_junction_filter", {})
    filter_by_junction = bool(junc_cfg.get("enabled", False))
    min_gc_at_junction = int(junc_cfg.get("min_gc_at_junction", 1))
    max_terminal_at_run = int(junc_cfg.get("max_terminal_at_run", 2))
    check_nick_dimer = bool(junc_cfg.get("check_nick_dimer", True))
    max_nick_dimer_bp = int(junc_cfg.get("max_nick_dimer_bp", 3))

    # Option B: Ligase Footprint Filter on BLAST Hits
    lig_footprint_cfg = arm_cfg.get("ligase_footprint_filter", {})
    enable_ligase_footprint = bool(lig_footprint_cfg.get("enabled", False))
    rescue_3prime_mismatches = bool(lig_footprint_cfg.get("rescue_3prime_mismatches", False))

    min_tm = float(arm_cfg.get("min_tm", arm_cfg.get("minArmTm", 40.0)))
    max_homopolymer = int(arm_cfg.get("maxHomopolymer", 4))
    max_hairpin_tm = float(arm_cfg.get("max_hairpin_tm", 30.0))
    max_arm_dimer_dg = float(arm_cfg.get("max_arm_dimer_dg", -1.0))

    blastdb_dir = arm_cfg.get("blastdbDir", "blastdb")
    threads = int(arm_cfg.get("threads", os.cpu_count() or 32))
    blast_min_ident = float(arm_cfg.get("blastMinIdentity", 70.0))
    max_allowed_hits = arm_cfg.get("maxAllowedHits", 0)

    print("=" * 70)
    print("STEP 4: PLP ARM GENERATION & DECOUPLED SPECIFICITY/THERMODYNAMIC SCREENING")
    print("=" * 70)
    print(f"Input Conserved FASTA  : {input_fasta}")
    print(f"Output CSV Target      : {output_csv}")
    print(f"Arm Sizes to Test      : {arm_lengths} nt (Footprint: {[2 * x for x in arm_lengths]} nt)")
    print(f"Threads                : {threads}")
    print("-" * 70)
    print("ACTIVE FILTERS:")
    print(f"  [ {'ON' if filter_by_tm else 'OFF'} ] Tm Filter         : >= {min_tm} °C (Primer3 default)")
    print(f"  [ {'ON' if filter_by_homopolymer else 'OFF'} ] Homopolymer Filter: Max Run < {max_homopolymer} nt")
    print(f"  [ {'ON' if filter_by_hairpin else 'OFF'} ] Hairpin Filter    : Max Hairpin Tm <= {max_hairpin_tm} °C (Primer3 default)")
    print(f"  [ {'ON' if filter_by_arm_dimer else 'OFF'} ] Arm Dimer Filter : ViennaRNA dG >= {max_arm_dimer_dg} kcal/mol")
    print(f"  [ {'ON' if filter_by_junction else 'OFF'} ] Junction Filter  : Min GC >= {min_gc_at_junction}, Max Term AT <= {max_terminal_at_run}")
    print(f"  [ {'ON' if enable_ligase_footprint else 'OFF'} ] Ligase Footprint: 3' Mismatch Rescue {'ENABLED' if rescue_3prime_mismatches else 'DISABLED'}")
    print("-" * 70)

    blast_dbs = discover_blast_dbs(blastdb_dir)

    if not os.path.exists(input_fasta):
        print(f"[Error] Input FASTA '{input_fasta}' not found.")
        return False

    records = list(SeqIO.parse(input_fasta, "fasta"))
    if not records:
        print(f"[Error] No sequences found in '{input_fasta}'.")
        return False

    # Identify on-target database if present in blastdb to avoid penalizing on-target binding
    target_db_name = None
    target_org = str(arm_cfg.get("target_organism") or cfg.get("target_organism") or "").strip().lower()
    if target_org:
        for db in blast_dbs:
            if db.lower() == target_org:
                target_db_name = db
                break
    if not target_db_name:
        search_ctx = (input_fasta + " " + (records[0].description if records else "")).lower()
        for db in blast_dbs:
            if db.lower() in search_ctx:
                target_db_name = db
                break

    off_target_dbs = [d for d in blast_dbs if d != target_db_name]

    print(f"Min Identity %         : >= {blast_min_ident}%")
    if isinstance(max_allowed_hits, (int, float)):
        print(f"Max Allowed Hits Rule  : Total off-target sum across background DBs <= {int(max_allowed_hits)}")
    else:
        print(f"Max Allowed Hits Rule  : {max_allowed_hits}")

    db_display = []
    if target_db_name:
        db_display.append(f"{target_db_name} (target)")
    db_display.extend(off_target_dbs)
    print(f"Databases Found        : {', '.join(db_display) if blast_dbs else 'None'}")
    print("=" * 70)

    candidates = []
    total_scanned = 0
    dropped_tm = 0
    dropped_homopolymer = 0
    dropped_hairpin = 0
    dropped_dimer = 0
    dropped_junction = 0

    print("\n[4.1] Scanning sliding windows & evaluating ViennaRNA thermodynamics...")

    for record in records:
        seq_str = str(record.seq).upper()
        target_len = len(seq_str)
        aln_block_start, _ = parse_aln_pos_offset(record.description)
        cons_meta = parse_conservation_metadata(record.description)

        for arm_len in arm_lengths:
            span = 2 * arm_len
            if target_len < span:
                continue

            for i in range(target_len - span + 1):
                total_scanned += 1
                local_start = i + 1
                local_end = i + span

                aln_start = aln_block_start + i
                aln_end = aln_start + span - 1
                aln_pos = f"{aln_start}-{aln_end}"

                target_region = seq_str[i : i + span]
                if any(base not in "ACGT" for base in target_region):
                    continue

                left_target = target_region[:arm_len]
                right_target = target_region[arm_len:]

                left_rc = str(Seq(left_target).reverse_complement())
                right_rc = str(Seq(right_target).reverse_complement())

                # 1. Tm Filter
                left_gc, left_tm = compute_metrics(left_rc)
                right_gc, right_tm = compute_metrics(right_rc)

                if filter_by_tm:
                    if left_tm < min_tm or right_tm < min_tm:
                        dropped_tm += 1
                        continue

                # 2. Homopolymer Filter
                left_homo_len, left_homo_seq = get_max_homopolymer(left_rc)
                right_homo_len, right_homo_seq = get_max_homopolymer(right_rc)

                if filter_by_homopolymer:
                    if left_homo_len >= max_homopolymer or right_homo_len >= max_homopolymer:
                        dropped_homopolymer += 1
                        continue

                # 3. Hairpin Filter
                left_hp_found, left_hp_tm, _ = check_hairpin(left_rc)
                right_hp_found, right_hp_tm, _ = check_hairpin(right_rc)

                if filter_by_hairpin:
                    if (left_hp_found == "Yes" and left_hp_tm > max_hairpin_tm) or \
                       (right_hp_found == "Yes" and right_hp_tm > max_hairpin_tm):
                        dropped_hairpin += 1
                        continue

                # 4. ViennaRNA Arm Dimer Filter
                dimer_db, vienna_dg = eval_arm_dimer_vienna(left_rc, right_rc)

                if filter_by_arm_dimer:
                    if vienna_dg < max_arm_dimer_dg:
                        dropped_dimer += 1
                        continue

                # 5. Ligation Junction Evaluation (Option A)
                junc_ok, junc_reason, junc_metrics = evaluate_ligation_junction(
                    left_rc,
                    right_rc,
                    min_gc_at_junction=min_gc_at_junction,
                    max_terminal_at_run=max_terminal_at_run,
                    check_nick_dimer=check_nick_dimer,
                    max_nick_dimer_bp=max_nick_dimer_bp,
                )
                if filter_by_junction and not junc_ok:
                    dropped_junction += 1
                    continue

                candidates.append({
                    "target_id": record.id,
                    "target_start": local_start,
                    "target_end": local_end,
                    "aln_pos": aln_pos,
                    "arm_len": arm_len,
                    "target_region": target_region,
                    "left_target": left_target,
                    "right_target": right_target,
                    "left_rc": left_rc,
                    "right_rc": right_rc,
                    "left_gc": left_gc,
                    "left_tm": left_tm,
                    "right_gc": right_gc,
                    "right_tm": right_tm,
                    "left_homo": f"{left_homo_len} ({left_homo_seq})",
                    "right_homo": f"{right_homo_len} ({right_homo_seq})",
                    "left_hp_tm": left_hp_tm,
                    "right_hp_tm": right_hp_tm,
                    "dimer_vienna_dg": vienna_dg,
                    "dimer_mfold_dg": "N/A",
                    "dimer_db": dimer_db,
                    "junc_gc": junc_metrics.get("Junction_GC", 0),
                    "junc_3p": junc_metrics.get("Junction_3prime_Base", ""),
                    "junc_5p": junc_metrics.get("Junction_5prime_Base", ""),
                    "conservation_pct": cons_meta.get("conservation_pct"),
                    "avg_conservation_pct": cons_meta.get("avg_conservation_pct"),
                    "organism": cons_meta.get("organism", ""),
                })

    print("-" * 70)
    print("THERMODYNAMIC & STRUCTURAL FILTERING BREAKDOWN:")
    print(f"  Total Windows Evaluated           : {total_scanned:,}")
    if filter_by_tm:
        pct_tm = (dropped_tm / total_scanned * 100) if total_scanned else 0
        print(f"  [-] Dropped by Tm (< {min_tm} °C)        : {dropped_tm:,} ({pct_tm:.2f}%)")
    if filter_by_homopolymer:
        pct_homo = (dropped_homopolymer / total_scanned * 100) if total_scanned else 0
        print(f"  [-] Dropped by Homopolymer (>= {max_homopolymer}) : {dropped_homopolymer:,} ({pct_homo:.2f}%)")
    if filter_by_hairpin:
        pct_hp = (dropped_hairpin / total_scanned * 100) if total_scanned else 0
        print(f"  [-] Dropped by Hairpin (> {max_hairpin_tm} °C)     : {dropped_hairpin:,} ({pct_hp:.2f}%)")
    if filter_by_arm_dimer:
        pct_dimer = (dropped_dimer / total_scanned * 100) if total_scanned else 0
        print(f"  [-] Dropped by Arm Dimer (< {max_arm_dimer_dg} kcal): {dropped_dimer:,} ({pct_dimer:.2f}%)")
    if filter_by_junction:
        pct_junc = (dropped_junction / total_scanned * 100) if total_scanned else 0
        print(f"  [-] Dropped by Junction Filter    : {dropped_junction:,} ({pct_junc:.2f}%)")

    surviving_pct = (len(candidates) / total_scanned * 100) if total_scanned else 0
    print(f"  [=] Surviving Thermodynamic Gate  : {len(candidates):,} ({surviving_pct:.2f}%)")
    print("-" * 70)

    if not candidates:
        print("[Warning] No candidates passed thermodynamic filters.")
        return False

    print("\n[4.2] Running BLAST specificity checks on contiguous footprint...")
    db_drop_counts = {db: 0 for db in blast_dbs}
    dropped_blast_total = 0

    if blast_dbs:
        unique_targets = {cand["target_region"]: f"T_{idx}" for idx, cand in enumerate(candidates)}
        db_hits = {db: {} for db in blast_dbs}
        db_lig_hits = {db: {} for db in blast_dbs} if enable_ligase_footprint else None
        footprint_len = 2 * arm_lengths[0]

        for db_name, db_prefix in blast_dbs.items():
            print(f"  -> Querying database: {db_name}...")
            if enable_ligase_footprint:
                hits, lig_hits = run_blastn_short_parallel(
                    unique_targets,
                    db_prefix,
                    arm_len=footprint_len,
                    min_identity=blast_min_ident,
                    threads=threads,
                    check_ligation_junction=True,
                )
                db_hits[db_name] = hits
                db_lig_hits[db_name] = lig_hits
            else:
                hits = run_blastn_short_parallel(
                    unique_targets,
                    db_prefix,
                    arm_len=footprint_len,
                    min_identity=blast_min_ident,
                    threads=threads,
                    check_ligation_junction=False,
                )
                db_hits[db_name] = hits

        # Auto-detect on-target organism if present in blastdb and not explicitly set
        if not target_db_name and blast_dbs and candidates:
            for db_name in blast_dbs:
                match_count = sum(1 for c in candidates if db_hits[db_name].get(c["target_region"], 0) > 0)
                if match_count / len(candidates) >= 0.70:
                    target_db_name = db_name
                    off_target_dbs = [d for d in blast_dbs if d != target_db_name]
                    print(f"  [On-Target Identified]: '{target_db_name}' matches {match_count}/{len(candidates)} ({match_count/len(candidates)*100:.1f}%) candidate windows.")
                    print(f"  [Off-Target Panel]   : {', '.join(off_target_dbs) if off_target_dbs else 'None'}")
                    break

        filtered_candidates = []
        rescued_blast_total = 0
        for cand in candidates:
            target_seq = cand["target_region"]
            passes_blast = True
            extra_counts = []
            lig_counts_list = []

            for db_name in blast_dbs:
                total_hits = db_hits[db_name].get(target_seq, 0)
                extra_counts.append(total_hits)
                if enable_ligase_footprint and db_lig_hits:
                    lig_counts_list.append(db_lig_hits[db_name].get(target_seq, 0))

            if isinstance(max_allowed_hits, (int, float)):
                # Sum of hits across all off-target organisms in blastdb
                off_target_sum = sum(
                    db_hits[db_name].get(target_seq, 0)
                    for db_name in off_target_dbs
                )
                if off_target_sum > int(max_allowed_hits):
                    # Check if 3' junction mismatch rescue is enabled
                    if enable_ligase_footprint and rescue_3prime_mismatches and db_lig_hits:
                        lig_off_target_sum = sum(
                            db_lig_hits[db_name].get(target_seq, 0)
                            for db_name in off_target_dbs
                        )
                        if lig_off_target_sum <= int(max_allowed_hits):
                            passes_blast = True
                            rescued_blast_total += 1
                            cand["rescued_by_3prime_mismatch"] = True
                        else:
                            passes_blast = False
                            dropped_blast_total += 1
                            for db_name in off_target_dbs:
                                if db_hits[db_name].get(target_seq, 0) > 0:
                                    db_drop_counts[db_name] += 1
                    else:
                        passes_blast = False
                        dropped_blast_total += 1
                        for db_name in off_target_dbs:
                            if db_hits[db_name].get(target_seq, 0) > 0:
                                db_drop_counts[db_name] += 1
            elif isinstance(max_allowed_hits, dict):
                for db_name in blast_dbs:
                    total_hits = db_hits[db_name].get(target_seq, 0)
                    allowed = max_allowed_hits.get(db_name)
                    if allowed is not None and total_hits > allowed:
                        if enable_ligase_footprint and rescue_3prime_mismatches and db_lig_hits:
                            lig_hits = db_lig_hits[db_name].get(target_seq, 0)
                            if lig_hits <= allowed:
                                continue
                        passes_blast = False
                        db_drop_counts[db_name] += 1
                if not passes_blast:
                    dropped_blast_total += 1

            if passes_blast:
                cand["blast_counts"] = extra_counts
                if enable_ligase_footprint:
                    cand["lig_counts"] = lig_counts_list
                filtered_candidates.append(cand)

        print("-" * 70)
        print("BLAST SPECIFICITY FILTERING BREAKDOWN:")
        print(f"  Total Candidates Submitted to BLAST: {len(candidates):,}")
        for db_name in blast_dbs:
            if db_name == target_db_name:
                continue
            drops = db_drop_counts.get(db_name, 0)
            print(f"  [-] Violations in {db_name:<16}: {drops:,}")

        if enable_ligase_footprint and rescue_3prime_mismatches:
            print(f"  [+] Rescued by 3' Nick Mismatch    : {rescued_blast_total:,}")

        blast_drop_pct = (dropped_blast_total / len(candidates) * 100) if candidates else 0
        surv_blast_pct = (len(filtered_candidates) / len(candidates) * 100) if candidates else 0
        print(f"  [-] Total Dropped by Specificity    : {dropped_blast_total:,} ({blast_drop_pct:.2f}%)")
        print(f"  [=] Confirmed Specific Candidates   : {len(filtered_candidates):,} ({surv_blast_pct:.2f}%)")
        print("-" * 70)

        candidates = filtered_candidates
    else:
        for cand in candidates:
            cand["blast_counts"] = []
            if enable_ligase_footprint:
                cand["lig_counts"] = []

    # -----------------------------------------------------------------
    # [4.3] Post-BLAST DNA mfold Profiling on Surviving Specific Candidates
    # -----------------------------------------------------------------
    print(f"\n[4.3] Profiling DNA secondary structure and mfold dG on {len(candidates):,} confirmed candidates...")
    for cand in candidates:
        mfold_dg = eval_arm_dimer_mfold(cand["left_rc"], cand["right_rc"])
        cand["dimer_mfold_dg"] = mfold_dg if mfold_dg is not None else "N/A"
    print(f"  [=] Completed mfold profiling for {len(candidates):,} candidates.")

    # -----------------------------------------------------------------
    # [4.4] CSV Export
    # -----------------------------------------------------------------
    print(f"\n[4.4] Writing qualified candidates to '{output_csv}'...")

    headers = [
        "Candidate_ID",
        "Target_ID",
        "aln_pos",
        "Arm_Length",
        "Target_Footprint",
        "Left_Target",
        "Right_Target",
        "Left_Arm",
        "Right_Arm",
        "Left_Arm_GC%",
        "Right_Arm_GC%",
        "Left_Arm_Tm",
        "Right_Arm_Tm",
        "Left_Arm_Homopolymer",
        "Right_Arm_Homopolymer",
        "Left_Arm_Hairpin_Tm",
        "Right_Arm_Hairpin_Tm",
        "Arm_Dimer_Vienna_dG",
        "Arm_Dimer_mfold_dG",
        "Junction_GC",
        "Junction_3prime_Base",
        "Junction_5prime_Base",
        "Conservation_%",
        "Avg_Conservation_%",
    ]
    for db_name in blast_dbs:
        headers.append(f"Target_Hits_{db_name}")
    if enable_ligase_footprint:
        for db_name in blast_dbs:
            headers.append(f"Ligation_Competent_Hits_{db_name}")

    out_dir = os.path.dirname(os.path.abspath(output_csv))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(output_csv, "w", newline="") as f_out:
        writer = csv.writer(f_out)
        writer.writerow(headers)
        for idx, cand in enumerate(candidates, start=1):
            row = [
                f"Candidate_{idx}",
                cand["target_id"],
                cand["aln_pos"],
                cand["arm_len"],
                cand["target_region"],
                cand["left_target"],
                cand["right_target"],
                cand["left_rc"],
                cand["right_rc"],
                cand["left_gc"],
                cand["right_gc"],
                cand["left_tm"],
                cand["right_tm"],
                cand["left_homo"],
                cand["right_homo"],
                cand["left_hp_tm"],
                cand["right_hp_tm"],
                cand["dimer_vienna_dg"],
                cand["dimer_mfold_dg"],
                cand.get("junc_gc", "N/A"),
                cand.get("junc_3p", "N/A"),
                cand.get("junc_5p", "N/A"),
                cand.get("conservation_pct") if cand.get("conservation_pct") is not None else "N/A",
                cand.get("avg_conservation_pct") if cand.get("avg_conservation_pct") is not None else "N/A",
            ] + cand["blast_counts"] + cand.get("lig_counts", [])
            writer.writerow(row)

    print(f"[Success] Exported {len(candidates):,} candidate arm pairs to {output_csv}\n")
    return True
