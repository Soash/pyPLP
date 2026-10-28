import csv
import multiprocessing as mp
import os
import sys
import time

from pyPLP.blast.discovery import discover_blast_dbs
from pyPLP.blast.search import run_blastn_short_parallel
from pyPLP.plp_generation.arm_sampler import sample_candidate_arms
from pyPLP.plp_generation.vienna_plot import render_vienna_svg
from pyPLP.plp_generation.worker import init_worker, worker_screen_batch
from pyPLP.thermodynamics.mfold_detailed import eval_structure_mfold_detailed


def generate_plp_for_arm(
    left_arm: str,
    right_arm: str,
    arm_meta: dict,
    cfg: dict,
    structures_dir: str,
    blast_dbs: dict,
    pool: mp.Pool,
    total_counter: mp.Value,
    lock: mp.Lock,
    stop_event: mp.Event | None = None,
) -> list[dict]:
    """Generates confirmed open, specific PLP probes for a single arm pair."""
    arm_cfg = cfg.get("arm_design", {})
    target_count = int(cfg.get("plp_per_arm", cfg.get("target_plp_per_arm", cfg.get("target_plp_count", 5))))
    threads = int(cfg.get("threads", arm_cfg.get("threads", os.cpu_count() or 32)))
    max_homo = int(cfg.get("maxHomopolymer", arm_cfg.get("maxHomopolymer", 4)))
    max_hairpin_tm = float(cfg.get("max_hairpin_tm", arm_cfg.get("max_hairpin_tm", 35.0)))
    min_plp_mfe_dg = float(cfg.get("min_plp_mfe_dg", -6.0))
    max_paired_arm_bases = int(cfg.get("max_paired_arm_bases", 4))
    max_paired_total_bases = int(cfg.get("max_paired_total_bases", 16))
    blast_min_ident = float(cfg.get("blastMinIdentity", arm_cfg.get("blastMinIdentity", 70.0)))
    max_allowed_hits = cfg.get("maxAllowedHits", arm_cfg.get("maxAllowedHits", 0))
    keep_ct_ps = bool(cfg.get("keep_ct_ps", False))

    len_l = len(left_arm)
    len_r = len(right_arm)
    avg_arm_len = (len_l + len_r) // 2

    arch_cfg = cfg.get("backbone_architecture") or {}
    spacer1_len = int(arch_cfg.get("spacer1_len", 4))
    rca_primer_binding_len = int(arch_cfg.get("rca_primer_binding_len", 20))
    spacer2_len = int(arch_cfg.get("spacer2_len", 4))
    probe_binding_len = int(arch_cfg.get("probe_binding_len", 20))
    spacer3_len = int(arch_cfg.get("spacer3_len", 0))
    rest_len = int(arch_cfg.get("rest_len", 4))
    re_site_str = str(arch_cfg.get("restriction_site", "")).strip()
    from pyPLP.plp_generation.segmented_backbone import resolve_restriction_site
    re_seq = resolve_restriction_site(re_site_str) if re_site_str else ""
    re_len = len(re_seq)
    min_tm_default = float(arm_cfg.get("min_tm", 40.0))
    min_rca_tm = float(arch_cfg.get("min_rca_primer_tm", min_tm_default))
    min_probe_tm = float(arch_cfg.get("min_probe_tm", min_tm_default))

    arch_bb_len = spacer1_len + rca_primer_binding_len + spacer2_len + probe_binding_len + spacer3_len + re_len + rest_len
    avg_arm = (len(left_arm) + len(right_arm)) // 2

    bb_factor = cfg.get("backbone_length_factor")
    if "backbone_length" in cfg and cfg["backbone_length"]:
        target_bb_len = int(cfg["backbone_length"])
    elif bb_factor is not None and str(bb_factor).strip():
        target_bb_len = int(bb_factor) * avg_arm
    else:
        target_bb_len = arch_bb_len

    if target_bb_len < arch_bb_len:
        target_bb_len = arch_bb_len
    elif target_bb_len > arch_bb_len:
        rest_len += (target_bb_len - arch_bb_len)
        arch_bb_len = target_bb_len

    backbone_len = target_bb_len

    effective_arch_cfg = {
        "spacer1_len": spacer1_len,
        "rca_primer_binding_len": rca_primer_binding_len,
        "spacer2_len": spacer2_len,
        "probe_binding_len": probe_binding_len,
        "spacer3_len": spacer3_len,
        "rest_len": rest_len,
        "restriction_site": re_site_str,
        "re_mode": arch_cfg.get("re_mode", "standalone"),
        "re_probe_offset": int(arch_cfg.get("re_probe_offset", 2)),
        "universal_rca_primer_seq": str(arch_cfg.get("universal_rca_primer_seq", "")).strip(),
        "universal_probe_seq": str(arch_cfg.get("universal_probe_seq", "")).strip(),
        "min_rca_primer_tm": min_rca_tm,
        "min_probe_tm": min_probe_tm,
        "screen_against": cfg.get("screen_against", []),
        "max_cross_dimer_dg": cfg.get("max_cross_dimer_dg", -5.0),
        "crispr_mode": bool(arch_cfg.get("crispr_mode", False) or cfg.get("crispr_mode", False)),
        "rest_tail": str(arch_cfg.get("rest_tail", "") or cfg.get("rest_tail", "")),
    }

    should_plot = bool(cfg.get("generate_plots", True)) and bool(structures_dir)
    if should_plot:
        vienna_dir = os.path.join(structures_dir, "vienna")
        mfold_dir = os.path.join(structures_dir, "mfold")
        os.makedirs(vienna_dir, exist_ok=True)
        os.makedirs(mfold_dir, exist_ok=True)
    else:
        vienna_dir = None
        mfold_dir = None

    worker_args = (
        left_arm,
        right_arm,
        backbone_len,
        max_homo,
        max_hairpin_tm,
        min_plp_mfe_dg,
        max_paired_arm_bases,
        max_paired_total_bases,
        effective_arch_cfg,
    )

    def _check_cand_blast(cand_to_test):
        if not blast_dbs:
            cand_to_test["blast_hits"] = []
            return True
        bb_seq = cand_to_test["bb"]
        q_dict = {bb_seq: "cand_q"}
        e_hits = []
        for db_n, db_p in blast_dbs.items():
            res_hits = run_blastn_short_parallel(
                q_dict,
                db_p,
                arm_len=backbone_len,
                min_identity=blast_min_ident,
                threads=min(4, threads),
            )
            e_hits.append(res_hits.get(bb_seq, 0))

        if isinstance(max_allowed_hits, (int, float)):
            if sum(e_hits) > int(max_allowed_hits):
                return False
        elif isinstance(max_allowed_hits, dict):
            for db_n, h_cnt in zip(blast_dbs, e_hits):
                allow = max_allowed_hits.get(db_n)
                if allow is not None and h_cnt > allow:
                    return False
        cand_to_test["blast_hits"] = e_hits
        return True

    valid_records = []
    arm_id = arm_meta.get("Candidate_ID", "Arm")

    if stop_event is not None:
        stop_event.clear()

    # Pre-populate parallel workers across all available threads
    pending_tasks = [pool.apply_async(worker_screen_batch, (worker_args,)) for _ in range(threads)]

    while len(valid_records) < target_count:
        new_pending = []
        for task in pending_tasks:
            if task.ready():
                try:
                    cand = task.get()
                except Exception:
                    cand = None

                if cand is not None:
                    if _check_cand_blast(cand):
                        valid_records.append(cand)
                        sys.stdout.write(
                            f"\r>> [{arm_id}] Optimizations Processed: {total_counter.value:,} | Confirmed Open: {len(valid_records)}/{target_count}"
                        )
                        sys.stdout.flush()
                        if len(valid_records) >= target_count:
                            break
                # Replenish with a fresh task if still needed
                if len(valid_records) < target_count:
                    new_pending.append(pool.apply_async(worker_screen_batch, (worker_args,)))
            else:
                new_pending.append(task)

        pending_tasks = new_pending

        if len(valid_records) >= target_count:
            break

        sys.stdout.write(
            f"\r>> [{arm_id}] Optimizations Processed: {total_counter.value:,} | Confirmed Open: {len(valid_records)}/{target_count}"
        )
        sys.stdout.flush()
        time.sleep(0.05)

    if stop_event is not None:
        stop_event.set()

    sys.stdout.write(
        f"\r>> [{arm_id}] Optimizations Processed: {total_counter.value:,} | Confirmed Open: {len(valid_records)}/{target_count}\n"
    )
    sys.stdout.flush()

    # Generate visual plots and detailed mfold structures for confirmed probes
    for idx, rec in enumerate(valid_records, start=1):
        plp_id = f"{arm_id}_PLP_{idx}"
        rec["plp_id"] = plp_id
        rec["left_arm"] = left_arm
        rec["right_arm"] = right_arm

        if should_plot and vienna_dir:
            svg_file = render_vienna_svg(
                plp_id, rec["full_plp"], rec["dot_bracket"], vienna_dir
            )
            rec["vienna_plot"] = svg_file

            m_dg, m_db, m_plot = eval_structure_mfold_detailed(
                rec["full_plp"],
                seq_id=plp_id,
                output_dir=mfold_dir,
                temperature=37.0,
                keep_ct_ps=keep_ct_ps,
            )
            rec["mfold_dg"] = m_dg if m_dg is not None else 0.0
            rec["mfold_db"] = m_db
            rec["mfold_plot"] = m_plot
        else:
            rec["vienna_plot"] = ""
            rec["mfold_dg"] = 0.0
            rec["mfold_db"] = ""
            rec["mfold_plot"] = ""

    return valid_records


def generate_plps_from_candidates(
    candidate_csv: str,
    output_csv: str,
    cfg: dict,
    job_dir: str = ".",
) -> bool:
    """Step 5: Samples N candidate arm pairs from candidate_plp_arms.csv and designs

    high-openness padlock probes with optimized, non-cross-reactive backbones.
    """
    raw_gen = cfg.get("plp_generation", cfg)
    if not raw_gen.get("enabled", True):
        print("[Notice] Step 5 (PLP Backbone Generation) is disabled in config.")
        return True

    arm_cfg = cfg.get("arm_design", {})
    gen_cfg = dict(raw_gen)
    for k in [
        "max_hairpin_tm",
        "maxHomopolymer",
        "blastdbDir",
        "threads",
        "blastMinIdentity",
        "maxAllowedHits",
    ]:
        if k not in gen_cfg:
            if k in arm_cfg:
                gen_cfg[k] = arm_cfg[k]
            elif k in cfg:
                gen_cfg[k] = cfg[k]

    sample_count = int(gen_cfg.get("sample_arms_count", 5))
    target_per_arm = int(gen_cfg.get("plp_per_arm", gen_cfg.get("target_plp_per_arm", 5)))
    threads = int(gen_cfg.get("threads", os.cpu_count() or 32))
    structures_dir = os.path.join(job_dir, gen_cfg.get("structures_dir", "plp_structures"))
    blastdb_dir = gen_cfg.get("blastdbDir", "blastdb")

    print("=" * 75)
    print("STEP 5: FULL PADLOCK PROBE (PLP) BACKBONE DESIGN & HIGH-OPENNESS OPTIMIZATION")
    print("=" * 75)
    print(f"Input Candidate Arms CSV : {candidate_csv}")
    print(f"Output PLP Probes CSV    : {output_csv}")
    print(f"Sample Arms Count (N)    : {sample_count}")
    print(f"Screening Engine         : RNAfold (ViennaRNA Mathews 2004 DNA model)")
    print(f"Profiling Engine         : DNA mfold (post-generation profiling)")
    print(f"Threads Allocated        : {threads}")
    print(f"Structures Directory     : {structures_dir}")
    print("-" * 75)

    if not os.path.isfile(candidate_csv):
        print(f"[Error] Candidate arms file '{candidate_csv}' not found.")
        return False

    sampled_arms = sample_candidate_arms(candidate_csv, sample_count=sample_count)
    if not sampled_arms:
        print("[Warning] No candidate arms available for PLP generation.")
        return False

    print(f"Selected {len(sampled_arms)} candidate arm pair(s) for backbone assembly:")
    for arm in sampled_arms:
        print(f"  -> {arm.get('Candidate_ID')} (Target: {arm.get('Target_ID')}, Pos: {arm.get('aln_pos')}, Left: {arm.get('Left_Arm')}, Right: {arm.get('Right_Arm')})")
    print("-" * 75)

    blast_dbs = discover_blast_dbs(blastdb_dir)
    print(f"Active Background DBs    : {', '.join(blast_dbs.keys()) if blast_dbs else 'None'}")
    print("=" * 75 + "\n")

    total_counter = mp.Value("L", 0)
    lock = mp.Lock()

    headers = [
        "PLP_ID",
        "Candidate_ID",
        "Target_ID",
        "aln_pos",
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
        "Restriction_Site_Seq",
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
    for db_name in blast_dbs:
        headers.append(f"Backbone_Hits_{db_name}")

    all_plp_records = []
    stop_event = mp.Event()

    with mp.Pool(
        processes=threads,
        initializer=init_worker,
        initargs=(total_counter, lock, stop_event),
    ) as pool:
        for idx, arm_meta in enumerate(sampled_arms, start=1):
            left_arm = arm_meta.get("Left_Arm") or arm_meta.get("Probe_Left_Arm_RC", "")
            right_arm = arm_meta.get("Right_Arm") or arm_meta.get("Probe_Right_Arm_RC", "")
            arm_id = arm_meta.get("Candidate_ID", f"Candidate_{idx}")

            print(f"\n[{idx}/{len(sampled_arms)}] Designing {target_per_arm} open backbones for {arm_id}...")

            records = generate_plp_for_arm(
                left_arm=left_arm,
                right_arm=right_arm,
                arm_meta=arm_meta,
                cfg=gen_cfg,
                structures_dir=structures_dir,
                blast_dbs=blast_dbs,
                pool=pool,
                total_counter=total_counter,
                lock=lock,
                stop_event=stop_event,
            )

            for r in records:
                r["Candidate_ID"] = arm_id
                r["Target_ID"] = arm_meta.get("Target_ID", "N/A")
                r["aln_pos"] = arm_meta.get("aln_pos", "N/A")
                all_plp_records.append(r)

    print("\n" + "=" * 75)
    print(f"PLP ASSEMBLY COMPLETE: Generated {len(all_plp_records)} open probes across {len(sampled_arms)} arm candidates.")
    print("=" * 75)

    os.makedirs(os.path.dirname(os.path.abspath(output_csv)), exist_ok=True)
    with open(output_csv, "w", newline="") as f_out:
        writer = csv.writer(f_out)
        writer.writerow(headers)
        for rec in all_plp_records:
            row = [
                rec["plp_id"],
                rec["Candidate_ID"],
                rec["Target_ID"],
                rec["aln_pos"],
                rec["full_plp"],
                len(rec["full_plp"]),
                rec["full_gc"],
                rec["full_tm"],
                rec.get("left_arm", rec["full_plp"][:len(sampled_arms[0].get("Left_Arm", ""))]),
                rec["bb"],
                rec.get("right_arm", rec["full_plp"][-len(sampled_arms[0].get("Right_Arm", "")):]),
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
                rec.get("Restriction_Site_Seq", ""),
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
            ] + rec.get("blast_hits", [])
            writer.writerow(row)

    # Generate dedicated LFA order sheet CSV
    lfa_csv = os.path.join(os.path.dirname(os.path.abspath(output_csv)), "lfa_rca_oligos.csv")
    lfa_headers = [
        "Oligo_Name",
        "Sequence_5to3",
        "5_Modification",
        "Length_nt",
        "GC%",
        "Tm_C",
        "Restriction_Site",
        "Role",
        "Associated_PLP_ID",
    ]
    with open(lfa_csv, "w", newline="") as f_lfa:
        lfa_writer = csv.writer(f_lfa)
        lfa_writer.writerow(lfa_headers)
        for rec in all_plp_records:
            pid = rec["plp_id"]
            # 1. Intact Padlock Probe
            lfa_writer.writerow([
                pid,
                rec["full_plp"],
                "5'-Phos",
                len(rec["full_plp"]),
                rec["full_gc"],
                rec["full_tm"],
                rec.get("Restriction_Site_Seq") or "None",
                "Padlock Probe (PLP)",
                pid,
            ])
            # 2. 5'-Biotin RCA Primer
            if rec.get("RCA_Primer_Seq_5to3"):
                lfa_writer.writerow([
                    f"{pid}_RCA_Primer",
                    rec["RCA_Primer_Seq_5to3"],
                    "5'-Biotin",
                    len(rec["RCA_Primer_Seq_5to3"]),
                    rec.get("RCA_Primer_GC%", ""),
                    rec.get("RCA_Primer_Tm", ""),
                    "N/A",
                    "RCA Amplification Primer",
                    pid,
                ])
            # 3. 5'-FAM Detection Probe
            if rec.get("Detection_Probe_Seq_5to3"):
                lfa_writer.writerow([
                    f"{pid}_Detection_Probe",
                    rec["Detection_Probe_Seq_5to3"],
                    "5'-6-FAM",
                    len(rec["Detection_Probe_Seq_5to3"]),
                    rec.get("Detection_Probe_GC%", ""),
                    rec.get("Detection_Probe_Tm", ""),
                    "N/A",
                    "LFA Fluorescent Detection Probe",
                    pid,
                ])

    print(f"[Success] Exported {len(all_plp_records)} designed PLPs to '{output_csv}'")
    print(f"[Success] Exported LFA assay order sheet to '{lfa_csv}'")
    print(f"Vienna diagrams saved in           -> '{os.path.join(structures_dir, 'vienna')}/'")
    print(f"IDT mfold structures saved in      -> '{os.path.join(structures_dir, 'mfold')}/'")
    print("=" * 75 + "\n")
    return True
