import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import sys

from pyPLP.alignment.mafft import run_mafft
from pyPLP.alignment.strand_sense import resolve_target_strand_fasta
from pyPLP.arm_design.designer import design_plp_arms
from pyPLP.checkpoint.manager import CheckpointManager
from pyPLP.conservation.extract_blocks import extract_conserved_blocks
from pyPLP.multiplex.panel_optimizer import design_multiplex_panel
from pyPLP.plp_generation.generator import generate_plps_from_candidates


def _sanitize_pathogen_name(stem: str) -> str:
    name = stem
    for prefix in ("sequences_", "sequence_"):
        if name.lower().startswith(prefix):
            name = name[len(prefix):]
    return name.upper() if name else stem


def run_pipeline(
    job_name: str,
    fasta_file: str | list[str],
    config_file: str | None = None,
    force: bool = False,
    multiplex: int | None = None,
    min_conservation: float | None = None,
    target_strand: str | None = None,
    restriction_site: str | None = None,
    fixed_oligos: str | list | dict | None = None,
    mafft_threads: int | None = None,
    generate_plp: bool | None = None,
    blast_min_identity: float | None = None,
) -> bool:
    """Executes the full PLP generation pipeline (supports single or multi-pathogen inputs)."""
    unified_config = {}
    if config_file:
        cfg_path = Path(config_file)
        if cfg_path.is_file():
            with open(cfg_path, "r", encoding="utf-8") as f:
                unified_config = json.load(f)
        else:
            print(f"[Warning] Config file '{config_file}' not found. Using default parameters.")

    if fixed_oligos:
        unified_config.setdefault("multiplex", {})["fixed_oligos"] = fixed_oligos

    # Parse alignment parameters (fallback to root keys if not nested)
    aln_cfg = unified_config.get("alignment", unified_config)
    mafft_env = aln_cfg.get("mafft_env", "mafft_env")
    mafft_mode = aln_cfg.get("mafft_mode", "auto")

    # Resolve MAFFT threads: CLI argument takes precedence, then config, default 8
    if mafft_threads is not None:
        mafft_threads = int(mafft_threads)
    else:
        mafft_threads = int(aln_cfg.get("mafft_threads", 8))

    min_conservation_pct = (
        float(min_conservation)
        if min_conservation is not None
        else float(aln_cfg.get("min_conservation_pct", 100.0))
    )
    min_block_len_bp = int(aln_cfg.get("min_block_len_bp", 30))
    ignore_gaps = bool(aln_cfg.get("ignore_gaps", False))

    target_strand_opt = (
        str(target_strand).strip().lower()
        if target_strand is not None
        else str(aln_cfg.get("target_strand", "positive")).strip().lower()
    )

    # Resolve blastMinIdentity override
    if blast_min_identity is not None:
        unified_config.setdefault("arm_design", {})["blastMinIdentity"] = float(blast_min_identity)

    # Resolve generate_plp override
    if generate_plp is not None:
        unified_config.setdefault("plp_generation", {})["enabled"] = bool(generate_plp)

    # Resolve restriction site
    arch_cfg = unified_config.get("plp_generation", {}).get("backbone_architecture", {})
    re_site_opt = (
        str(restriction_site).strip()
        if restriction_site is not None
        else str(arch_cfg.get("restriction_site", "")).strip()
    )
    if restriction_site is not None:
        if "plp_generation" not in unified_config:
            unified_config["plp_generation"] = {}
        if "backbone_architecture" not in unified_config["plp_generation"]:
            unified_config["plp_generation"]["backbone_architecture"] = {}
        unified_config["plp_generation"]["backbone_architecture"]["restriction_site"] = re_site_opt

    job_path = Path(job_name)
    if not job_path.is_absolute():
        parts = job_path.parts
        if not parts or parts[0] != "results":
            job_path = Path("results") / job_path
    job_dir = job_path.resolve()
    input_dir = job_dir / "input"
    input_dir.mkdir(parents=True, exist_ok=True)

    ckpt = CheckpointManager(job_dir)

    if force:
        print("[Notice] --force enabled. Clearing existing checkpoints...")
        shutil.rmtree(ckpt.ckpt_dir, ignore_errors=True)
        ckpt.ckpt_dir.mkdir(parents=True, exist_ok=True)

    if isinstance(fasta_file, str):
        fasta_files = [fasta_file]
    else:
        fasta_files = list(fasta_file)

    plp_cfg = unified_config.get("plp_generation", {})
    arm_cfg = unified_config.get("arm_design", {})
    auto_fallback_setting = bool(aln_cfg.get("auto_fallback", True))
    min_fallback_setting = float(aln_cfg.get("min_fallback_pct", 65.0))
    fallback_step_setting = float(aln_cfg.get("fallback_step", 1.0))

    designed_plps_csv = job_dir / "designed_plp_probes.csv"
    multiplex_panel_csv = job_dir / "multiplex_panel_probes.csv"
    candidate_arms_csv = job_dir / "candidate_plp_arms.csv"

    print("=" * 80)
    print("pyPLP PIPELINE INITIALIZATION & ACTIVE FILTER CONFIGURATION")
    print("=" * 80)
    print(f"Job Name (Output Dir)  : {job_dir}")
    print(f"Input FASTA File(s)    : {', '.join(str(Path(f).name) for f in fasta_files)}")
    print(f"Config File Loaded     : {config_file if (config_file and os.path.isfile(config_file)) else ('Defaults (No config file provided)' if not config_file else config_file + ' [Not Found, using defaults]')}")
    print("-" * 80)
    print("[1. Multiple Sequence Alignment & Conservation Filters]")
    print(f"  - MAFFT CPU Threads    : {mafft_threads} (mode: {mafft_mode}, env: {mafft_env})")
    print(f"  - Genomic Strand Sense : {target_strand_opt}")
    print(f"  - Conservation Cutoff  : {min_conservation_pct:.1f}%")
    print(f"  - Minimum Block Length : {min_block_len_bp} bp")
    print(f"  - Automatic Fallback   : {'ENABLED (stepping down to ' + str(min_fallback_setting) + '%, step ' + str(fallback_step_setting) + '%)' if auto_fallback_setting else 'DISABLED (strictly requiring ' + str(min_conservation_pct) + '%)'}")
    print(f"  - Gap Handling Mode    : {'Ignore gaps' if ignore_gaps else 'Strict (gaps penalize conservation)'}")
    print("-" * 80)
    print("[2. PLP Arm Design & Specificity Filters]")
    print(f"  - Arm Lengths Tested   : {arm_cfg.get('armLengths', [15])} nt")
    print(f"  - Tm Filter            : [{'ON' if arm_cfg.get('filter_by_tm', True) else 'OFF'}] >= {arm_cfg.get('min_tm', 40.0)} °C")
    print(f"  - Homopolymer Filter   : [{'ON' if arm_cfg.get('filter_by_homopolymer', True) else 'OFF'}] Max Run < {arm_cfg.get('maxHomopolymer', 4)} nt")
    print(f"  - Hairpin Filter       : [{'ON' if arm_cfg.get('filter_by_hairpin', True) else 'OFF'}] Max Hairpin Tm <= {arm_cfg.get('max_hairpin_tm', 30.0)} °C")
    print(f"  - Arm Dimer Filter     : [{'ON' if arm_cfg.get('filter_by_arm_dimer', True) else 'OFF'}] ViennaRNA dG >= {arm_cfg.get('max_arm_dimer_dg', -1.0)} kcal/mol")
    print(f"  - BLAST Specificity    : [ON] Identity >= {arm_cfg.get('blastMinIdentity', 70.0):.1f}%, Max Allowed Hits = {arm_cfg.get('maxAllowedHits', 0)}")
    print(f"  - Ligation Junction    : [{'ON' if arm_cfg.get('ligation_junction_filter', {}).get('enabled', False) else 'OFF'}]")
    print(f"  - Ligase Footprint     : [{'ON' if arm_cfg.get('ligase_footprint_filter', {}).get('enabled', False) else 'OFF'}]")
    print("-" * 80)
    print("[3. Step 5 Backbone Assembly & Multiplexing]")
    print(f"  - Step 5 PLP Assembly  : {'ENABLED' if plp_cfg.get('enabled', False) else 'DISABLED (Terminating after Step 4 arm screening)'}")
    if multiplex or unified_config.get("multiplex", {}).get("target_channels"):
        print(f"  - Multiplex Panel Mode : ENABLED ({multiplex or unified_config.get('multiplex', {}).get('target_channels')}-plex panel)")
    print("=" * 80 + "\n")

    # =============================================================
    # Branch A: Multi-Pathogen Multiplex Pipeline (len(fasta_files) > 1)
    # =============================================================
    if len(fasta_files) > 1:
        print(f"\n[ pyPLP Multi-Pathogen Pipeline ] Processing {len(fasta_files)} pathogen targets...")
        all_candidates = []
        combined_rows = []
        pathogen_summary_list = []

        auto_fallback = bool(aln_cfg.get("auto_fallback", True))
        min_fallback_pct = float(aln_cfg.get("min_fallback_pct", 65.0))
        fallback_step = float(aln_cfg.get("fallback_step", 1.0))

        channel_targets = []
        plp_cfg = unified_config.get("plp_generation", {})
        sample_arms_count = int(plp_cfg.get("sample_arms_count", 1))

        for f_item in fasta_files:
            src_fpath = Path(f_item).resolve()
            if not src_fpath.is_file():
                print(f"[Error] Input FASTA not found: {src_fpath}", file=sys.stderr)
                return False

            raw_stem = src_fpath.stem
            pathogen_label = _sanitize_pathogen_name(raw_stem)
            print(f"\n--- Pathogen: {pathogen_label} ({src_fpath.name}) ---")

            dest_fasta = input_dir / f"{raw_stem}.fasta"
            aligned_fasta = job_dir / f"aligned_{raw_stem}.fasta"
            conserved_fasta = job_dir / f"conserved_blocks_{raw_stem}.fasta"
            cand_arms_csv = job_dir / f"candidate_plp_arms_{raw_stem}.csv"

            # Step 1: Input Setup
            step1 = f"step_1_input_{raw_stem}"
            if ckpt.is_completed(step1, dest_fasta):
                print(f"[Resume] Skipping {step1}: input already prepared.")
            else:
                if target_strand_opt == "positive":
                    shutil.copyfile(src_fpath, dest_fasta)
                else:
                    print(f"[{pathogen_label} Step 1] Preprocessing strand sense: '{target_strand_opt}'...")
                    resolve_target_strand_fasta(
                        str(src_fpath),
                        target_strand=target_strand_opt,
                        output_fasta=str(dest_fasta),
                    )
                ckpt.mark_completed(step1, {"input_file": str(dest_fasta), "target_strand": target_strand_opt})

            # Step 2: MAFFT Alignment
            step2 = f"step_2_mafft_{raw_stem}"
            if ckpt.is_completed(step2, aligned_fasta):
                print(f"[Resume] Skipping {step2}: alignment already exists.")
            else:
                print(f"[{pathogen_label} Step 2] Executing MAFFT alignment...")
                success = run_mafft(
                    input_fasta=dest_fasta,
                    output_fasta=aligned_fasta,
                    mafft_env=mafft_env,
                    mafft_mode=mafft_mode,
                    mafft_threads=mafft_threads,
                )
                if not success:
                    print(f"[Error] MAFFT alignment failed for {pathogen_label}.", file=sys.stderr)
                    return False
                ckpt.mark_completed(step2, {"output_file": str(aligned_fasta)})

            # Step 3 & 4: Conserved Extraction & Arm Screening with Coordinated Fallback
            step3 = f"step_3_conserved_{raw_stem}"
            step4 = f"step_4_arm_screening_{raw_stem}"

            if ckpt.is_completed(step4, cand_arms_csv) and ckpt.is_completed(step3, conserved_fasta):
                print(f"[Resume] Skipping {step3} and {step4}: candidate arms already designed.")
                meta = ckpt.get_metadata(step3)
                actual_pct = float(meta.get("actual_conservation_pct", min_conservation_pct))
                blocks_found = int(meta.get("blocks_found", 1))
            else:
                curr_target_pct = float(min_conservation_pct)
                actual_pct = curr_target_pct
                blocks_found = 0
                success_arms = False

                while curr_target_pct >= min_fallback_pct:
                    res = extract_conserved_blocks(
                        input_file=str(aligned_fasta),
                        output_file=str(conserved_fasta),
                        min_conservation_pct=curr_target_pct,
                        min_block_len_bp=min_block_len_bp,
                        ignore_gaps=ignore_gaps,
                        auto_fallback=False,
                        organism_name=pathogen_label,
                    )
                    if res and res.block_count > 0:
                        print(f"[{pathogen_label} Step 4] Designing arms on {res.block_count} conserved blocks at {curr_target_pct:.1f}%...")
                        arms_ok = design_plp_arms(
                            input_fasta=str(conserved_fasta),
                            output_csv=str(cand_arms_csv),
                            cfg=unified_config,
                        )
                        has_candidates = False
                        if arms_ok and os.path.isfile(cand_arms_csv):
                            with open(cand_arms_csv, "r", encoding="utf-8") as f:
                                has_candidates = sum(1 for line in f if line.strip()) > 1

                        if arms_ok and has_candidates:
                            actual_pct = curr_target_pct
                            blocks_found = res.block_count
                            success_arms = True
                            ckpt.mark_completed(step3, {
                                "output_file": str(conserved_fasta),
                                "actual_conservation_pct": actual_pct,
                                "blocks_found": blocks_found,
                            })
                            ckpt.mark_completed(step4, {"output_file": str(cand_arms_csv)})
                            print(f"[Success] Successfully designed qualified arms for {pathogen_label} at {actual_pct:.1f}% conservation.")
                            break
                        else:
                            if not auto_fallback:
                                break
                            print(f"[{pathogen_label}] Conserved blocks at {curr_target_pct:.1f}% did not yield qualified arms. Stepping down fallback...")
                    else:
                        if not auto_fallback:
                            break
                    curr_target_pct = round(curr_target_pct - fallback_step, 2)

                if not success_arms:
                    print(f"[Error] Arm design failed for {pathogen_label}: no qualified arms found down to {min_fallback_pct:.1f}%.", file=sys.stderr)
                    return False

            pathogen_summary_list.append({
                "Pathogen": pathogen_label,
                "FASTA_File": src_fpath.name,
                "Requested_Conservation_%": min_conservation_pct,
                "Actual_Conservation_%": actual_pct,
                "Conserved_Blocks_Found": blocks_found,
                "Fallback_Triggered": "Yes" if actual_pct < min_conservation_pct else "No",
            })

            # Read candidate arms for this pathogen
            pathogen_candidates = []
            if os.path.isfile(cand_arms_csv):
                with open(cand_arms_csv, "r", encoding="utf-8") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        row_copy = dict(row)
                        row_copy["Pathogen"] = pathogen_label
                        if not row_copy.get("Conservation_%") or row_copy.get("Conservation_%") == "N/A":
                            row_copy["Conservation_%"] = f"{actual_pct:.1f}"
                        combined_rows.append(row_copy)
                        pathogen_candidates.append({
                            "name": f"{pathogen_label}_{row.get('Candidate_ID', '')}",
                            "candidate_id": row.get("Candidate_ID", ""),
                            "target_id": row.get("Target_ID", ""),
                            "pathogen": pathogen_label,
                            "organism": pathogen_label,
                            "left_arm": row.get("Left_Arm"),
                            "right_arm": row.get("Right_Arm"),
                            "conservation_pct": row_copy.get("Conservation_%"),
                            "avg_conservation_pct": row_copy.get("Avg_Conservation_%"),
                        })

            if pathogen_candidates:
                all_candidates.append(pathogen_candidates[0])
                sampled_arms = pathogen_candidates[:sample_arms_count]
                channel_targets.append({
                    "pathogen": pathogen_label,
                    "arms": sampled_arms,
                })
            else:
                print(f"[Warning] No qualified candidate arms found for pathogen {pathogen_label}.", file=sys.stderr)

        # Write aggregated candidate_plp_arms.csv with Pathogen and Conservation_% up front
        if combined_rows:
            all_keys = set()
            for r in combined_rows:
                all_keys.update(r.keys())
            priority_cols = ["Pathogen", "Conservation_%", "Avg_Conservation_%"]
            fnames = [c for c in priority_cols if c in all_keys]
            for k in combined_rows[0].keys():
                if k not in fnames:
                    fnames.append(k)
            for k in all_keys:
                if k not in fnames:
                    fnames.append(k)
            for r in combined_rows:
                for k in fnames:
                    if k not in r:
                        r[k] = "N/A"
            with open(candidate_arms_csv, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=fnames)
                writer.writeheader()
                writer.writerows(combined_rows)

        # Write target_conservation_summary.csv
        summary_csv = job_dir / "target_conservation_summary.csv"
        if pathogen_summary_list:
            with open(summary_csv, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=list(pathogen_summary_list[0].keys()))
                writer.writeheader()
                writer.writerows(pathogen_summary_list)
            print(f"\n[Conservation] Saved organism conservation summary to: {summary_csv}")

        k = multiplex or unified_config.get("multiplex", {}).get("target_channels") or len(channel_targets)
        k = int(k)
        step5 = f"step_5_multiplex_multipathogen_{k}plex"
        if ckpt.is_completed(step5, multiplex_panel_csv):
            print(f"[Resume] Skipping {step5}: multiplex panel CSV already exists.")
        else:
            print(f"\n[Step 5] Triggering Orthogonal Multi-Pathogen ({k}-Plex) Panel Assembly...")
            if len(channel_targets) < 2:
                print(f"[Error] Need at least 2 target channels for multiplex panel, found {len(channel_targets)}.", file=sys.stderr)
                return False

            mplex_cfg = unified_config.get("multiplex", {})
            labels = mplex_cfg.get("channel_labels", ["5'-Biotin", "5'-DIG", "5'-TAMRA", "5'-Cy5", "5'-TexasRed"])
            univ_probe = mplex_cfg.get("universal_detection_probe", mplex_cfg.get("universal_probe", ""))
            cross_dg = float(mplex_cfg.get("max_cross_dimer_dg", -2.0))

            res = design_multiplex_panel(
                targets=channel_targets[:k],
                output_dir=str(job_dir),
                universal_detection_probe=univ_probe,
                channel_labels=labels[:k],
                max_cross_dimer_dg=cross_dg,
                restriction_site=re_site_opt,
                cfg=unified_config,
            )
            if not res:
                print("[Error] Step 5 Multi-pathogen multiplex panel generation failed.", file=sys.stderr)
                return False
            ckpt.mark_completed(step5, {"output_file": str(multiplex_panel_csv)})

        print(f"\n[Finished] Multi-Pathogen Job '{job_name}' completed successfully.")
        print(f"Aggregated Candidate Arms CSV: {candidate_arms_csv}")
        print(f"Multiplex Panel CSV          : {multiplex_panel_csv}")
        return True

    # =============================================================
    # Branch B: Single FASTA Pipeline (Preserving existing workflow)
    # =============================================================
    dest_fasta = input_dir / "sequences.fasta"
    aligned_fasta = job_dir / "aligned_sequences.fasta"
    out_conserved_name = f"conserved_blocks_{int(min_conservation_pct)}.fasta"
    conserved_fasta = job_dir / out_conserved_name

    # -------------------------------------------------------------
    # Step 1: Input Setup
    # -------------------------------------------------------------
    step1 = "step_1_input_setup"
    if ckpt.is_completed(step1, dest_fasta):
        print(f"[Resume] Skipping {step1}: input already prepared.")
    else:
        src_fasta = Path(fasta_files[0]).resolve()
        if target_strand_opt == "positive":
            print(f"[Step 1] Copying input FASTA to {dest_fasta}...")
            shutil.copyfile(src_fasta, dest_fasta)
        else:
            print(f"[Step 1] Preprocessing strand sense: '{target_strand_opt}' to {dest_fasta}...")
            resolve_target_strand_fasta(
                str(src_fasta),
                target_strand=target_strand_opt,
                output_fasta=str(dest_fasta),
            )
        ckpt.mark_completed(step1, {"input_file": str(dest_fasta), "target_strand": target_strand_opt})

    # -------------------------------------------------------------
    # Step 2: MAFFT Multiple Sequence Alignment
    # -------------------------------------------------------------
    step2 = "step_2_mafft_alignment"
    if ckpt.is_completed(step2, aligned_fasta):
        print(f"[Resume] Skipping {step2}: alignment already exists.")
    else:
        print(f"[Step 2] Executing MAFFT alignment...")
        success = run_mafft(
            input_fasta=dest_fasta,
            output_fasta=aligned_fasta,
            mafft_env=mafft_env,
            mafft_mode=mafft_mode,
            mafft_threads=mafft_threads,
        )
        if not success:
            print("[Error] MAFFT alignment step failed.", file=sys.stderr)
            return False
        ckpt.mark_completed(step2, {"output_file": str(aligned_fasta)})

    # -------------------------------------------------------------
    # Step 3: Extract Conserved Blocks
    # -------------------------------------------------------------
    auto_fallback = bool(aln_cfg.get("auto_fallback", True))
    min_fallback_pct = float(aln_cfg.get("min_fallback_pct", 65.0))
    fallback_step = float(aln_cfg.get("fallback_step", 1.0))

    step3 = "step_3_conserved_extraction"
    actual_pct = min_conservation_pct
    blocks_found = 0
    if ckpt.is_completed(step3, conserved_fasta):
        print(f"[Resume] Skipping {step3}: conserved blocks already exist.")
        meta = ckpt.get_metadata(step3)
        actual_pct = float(meta.get("actual_conservation_pct", min_conservation_pct))
        blocks_found = int(meta.get("blocks_found", 1))
    else:
        print(
            f"\n[Step 3] Extracting conserved blocks:"
            f"\n  -> Target Conservation : {min_conservation_pct:.1f}%"
            f"\n  -> Minimum Block Length: {min_block_len_bp} bp"
            f"\n  -> Auto-Fallback Search: {'ENABLED (stepping down to ' + str(min_fallback_pct) + '%, step ' + str(fallback_step) + '%)' if auto_fallback else 'DISABLED (strictly testing ' + str(min_conservation_pct) + '% only)'}"
            f"\n  -> Gap Handling Mode   : {'Ignore gaps' if ignore_gaps else 'Strict (gaps penalize conservation)'}"
        )
        res = extract_conserved_blocks(
            input_file=str(aligned_fasta),
            output_file=str(conserved_fasta),
            min_conservation_pct=min_conservation_pct,
            min_block_len_bp=min_block_len_bp,
            ignore_gaps=ignore_gaps,
            auto_fallback=auto_fallback,
            min_fallback_pct=min_fallback_pct,
            fallback_step=fallback_step,
            organism_name=Path(fasta_files[0]).stem,
        )
        if not res:
            print("[Error] No conserved blocks met the criteria. Pipeline stopped.", file=sys.stderr)
            return False
        actual_pct = res.actual_pct
        blocks_found = res.block_count
        ckpt.mark_completed(step3, {
            "output_file": str(conserved_fasta),
            "actual_conservation_pct": actual_pct,
            "blocks_found": blocks_found,
        })

    summary_csv = job_dir / "target_conservation_summary.csv"
    with open(summary_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "Pathogen", "FASTA_File", "Requested_Conservation_%", "Actual_Conservation_%", "Conserved_Blocks_Found", "Fallback_Triggered"
        ])
        writer.writeheader()
        writer.writerow({
            "Pathogen": Path(fasta_files[0]).stem,
            "FASTA_File": Path(fasta_files[0]).name,
            "Requested_Conservation_%": min_conservation_pct,
            "Actual_Conservation_%": actual_pct,
            "Conserved_Blocks_Found": blocks_found,
            "Fallback_Triggered": "Yes" if actual_pct < min_conservation_pct else "No",
        })
    print(f"\n[Conservation] Saved organism conservation summary to: {summary_csv}")

    # -------------------------------------------------------------
    # Step 4: Generate PLP Arms & Screening
    # -------------------------------------------------------------
    step4 = "step_4_plp_arm_screening"
    if ckpt.is_completed(step4, candidate_arms_csv):
        print(f"[Resume] Skipping {step4}: candidate PLP arms already designed.")
    else:
        print(f"\n[Step 4] Triggering PLP arm generation and BLAST screening...")
        success = design_plp_arms(
            input_fasta=str(conserved_fasta),
            output_csv=str(candidate_arms_csv),
            cfg=unified_config,
        )
        if not success:
            print("[Error] PLP arm design or specificity screening failed.", file=sys.stderr)
            return False
        ckpt.mark_completed(step4, {"output_file": str(candidate_arms_csv)})

    # -------------------------------------------------------------
    # Step 5: Full PLP Generation / Multiplex Panel Assembly
    # -------------------------------------------------------------
    num_mplex = multiplex or unified_config.get("multiplex", {}).get("target_channels")
    if num_mplex and int(num_mplex) > 1:
        k = int(num_mplex)
        step5 = f"step_5_multiplex_{k}plex_assembly"
        if ckpt.is_completed(step5, multiplex_panel_csv):
            print(f"[Resume] Skipping {step5}: multiplex panel CSV already exists.")
        else:
            print(f"\n[Step 5] Triggering Orthogonal Multiplex ({k}-Plex) Panel Assembly...")
            candidates = []
            seen_targets = set()
            with open(candidate_arms_csv, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    tgt_id = row.get("Target_ID", "Region")
                    if tgt_id not in seen_targets:
                        seen_targets.add(tgt_id)
                        candidates.append({
                            "name": f"{tgt_id}_{row.get('Candidate_ID', '')}",
                            "left_arm": row.get("Left_Arm"),
                            "right_arm": row.get("Right_Arm"),
                        })
                        if len(candidates) >= k:
                            break

            if len(candidates) < k:
                with open(candidate_arms_csv, "r", encoding="utf-8") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        cand_entry = {
                            "name": f"{row.get('Target_ID', 'R')}_{row.get('Candidate_ID', '')}",
                            "left_arm": row.get("Left_Arm"),
                            "right_arm": row.get("Right_Arm"),
                        }
                        if cand_entry not in candidates:
                            candidates.append(cand_entry)
                        if len(candidates) >= k:
                            break

            if len(candidates) < 2:
                print(f"[Error] Need at least 2 candidate arms for multiplex panel, found {len(candidates)}.", file=sys.stderr)
                return False

            mplex_cfg = unified_config.get("multiplex", {})
            labels = mplex_cfg.get("channel_labels", ["5'-Biotin", "5'-DIG", "5'-TAMRA"])
            univ_probe = mplex_cfg.get("universal_detection_probe", mplex_cfg.get("universal_probe", ""))
            cross_dg = float(mplex_cfg.get("max_cross_dimer_dg", -2.0))

            res = design_multiplex_panel(
                targets=candidates[:k],
                output_dir=str(job_dir),
                universal_detection_probe=univ_probe,
                channel_labels=labels[:k],
                max_cross_dimer_dg=cross_dg,
                restriction_site=re_site_opt,
                cfg=unified_config,
            )
            if not res:
                print("[Error] Step 5 Multiplex panel generation failed.", file=sys.stderr)
                return False
            ckpt.mark_completed(step5, {"output_file": str(multiplex_panel_csv)})
    elif plp_cfg.get("enabled", False):
        step5 = "step_5_plp_backbone_assembly"
        if ckpt.is_completed(step5, designed_plps_csv):
            print(f"[Resume] Skipping {step5}: designed PLPs CSV already exists.")
        else:
            print(f"\n[Step 5] Triggering full PLP assembly and backbone optimization...")
            success = generate_plps_from_candidates(
                candidate_csv=str(candidate_arms_csv),
                output_csv=str(designed_plps_csv),
                cfg=unified_config,
                job_dir=str(job_dir),
            )
            if not success:
                print("[Error] Step 5 PLP generation failed.", file=sys.stderr)
                return False
            ckpt.mark_completed(step5, {"output_file": str(designed_plps_csv)})

    print(f"\n[Finished] Job '{job_name}' completed all steps successfully.")
    print(f"Candidate Arms CSV : {candidate_arms_csv}")
    if os.path.isfile(multiplex_panel_csv):
        print(f"Multiplex Panel CSV: {multiplex_panel_csv}")
    elif os.path.isfile(designed_plps_csv):
        print(f"Designed PLPs CSV  : {designed_plps_csv}")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="pyPLP: Automated MAFFT alignment, conservation extraction, and PLP arm design pipeline."
    )
    parser.add_argument("-j", "--job-name", required=True, help="Job name (saved under results/<job_name>)")
    parser.add_argument(
        "-f",
        "--fasta",
        required=True,
        nargs="+",
        help="Path to input raw FASTA file(s). Pass multiple files for multi-pathogen multiplex panel design.",
    )
    parser.add_argument("-c", "--config", default=None, help="Path to unified JSON config (optional)")
    parser.add_argument(
        "-t",
        "--threads",
        "--mafft-threads",
        dest="mafft_threads",
        type=int,
        default=None,
        help="Number of CPU cores for MAFFT alignment (default: 8)",
    )
    parser.add_argument(
        "-b",
        "--blast-min-identity",
        type=float,
        default=None,
        help="BLAST minimum sequence identity percentage threshold to flag off-target hits (default: 70.0)",
    )
    parser.add_argument(
        "--generate-plp",
        action="store_true",
        default=None,
        help="Enable Step 5 full PLP backbone assembly (default: False, arm screening only)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Ignore checkpoints and force re-execution of all steps",
    )
    parser.add_argument(
        "-m",
        "--multiplex",
        type=int,
        default=None,
        help="Assemble an orthogonal multiplex panel across N distinct targets (e.g. -m 3)",
    )
    parser.add_argument(
        "-p",
        "--min-conservation",
        type=float,
        default=None,
        help="Initial minimum conservation percentage (e.g. -p 100.0). If not found, falls back step-by-step.",
    )
    parser.add_argument(
        "--target-strand",
        choices=["positive", "negative", "both"],
        default=None,
        help="Genomic RNA sense: 'positive' (default, +ssRNA/mRNA), 'negative' (-ssRNA vRNA, auto-generates antigenomic cRNA), or 'both'",
    )
    parser.add_argument(
        "--restriction-site",
        default=None,
        help="Optional restriction enzyme cleavage site (e.g. BstUI, EcoRV, HaeIII, or raw ACGT sequence like CGCG) inserted into PLP backbone",
    )
    parser.add_argument(
        "--fixed-oligos",
        default=None,
        help="Path to CSV file or string containing fixed reference oligos (PLP, RCA_Primer, PLP_Detection_Probe)",
    )

    args = parser.parse_args()
    fasta_arg = args.fasta if len(args.fasta) > 1 else args.fasta[0]
    success = run_pipeline(
        job_name=args.job_name,
        fasta_file=fasta_arg,
        config_file=args.config,
        force=args.force,
        multiplex=args.multiplex,
        min_conservation=args.min_conservation,
        target_strand=args.target_strand,
        restriction_site=args.restriction_site,
        fixed_oligos=args.fixed_oligos,
        mafft_threads=args.mafft_threads,
        generate_plp=args.generate_plp,
        blast_min_identity=args.blast_min_identity,
    )
    if not success:
        sys.exit(1)


if __name__ == "__main__":
    main()
