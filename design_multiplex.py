#!/usr/bin/env python3
"""design_multiplex.py: Standalone CLI and Python tool to design fully orthogonal,

non-cross-reactive multiplex (triplex / K-plex) padlock probe panels.
"""

import argparse
import csv
import json
import os
import sys
from pathlib import Path

# Ensure pyPLP is importable
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

from pyPLP.multiplex.panel_optimizer import design_multiplex_panel


def parse_targets_csv(csv_path: str) -> list[dict]:
    """Parses a CSV containing Target_Name (or Name), Left_Arm, Right_Arm."""
    targets = []
    if not os.path.isfile(csv_path):
        raise FileNotFoundError(f"Targets CSV file not found: '{csv_path}'")

    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for idx, row in enumerate(reader, start=1):
            name = row.get("Target_Name") or row.get("Target") or row.get("Name") or f"Target_{idx}"
            left = row.get("Left_Arm") or row.get("Left") or ""
            right = row.get("Right_Arm") or row.get("Right") or ""
            if not left or not right:
                continue
            targets.append({
                "name": name.strip(),
                "left_arm": left.strip().upper(),
                "right_arm": right.strip().upper(),
            })

    if not targets:
        raise ValueError(f"No valid targets found in '{csv_path}'. Required columns: Target_Name, Left_Arm, Right_Arm.")
    return targets


def parse_inline_targets(target_str: str) -> list[dict]:
    """Parses an inline comma-separated string formatted as: Name:LeftArm:RightArm,Name2:..."""
    targets = []
    tokens = [t.strip() for t in target_str.split(",") if t.strip()]
    for idx, tok in enumerate(tokens, start=1):
        parts = [p.strip() for p in tok.split(":") if p.strip()]
        if len(parts) == 3:
            name, left, right = parts
        elif len(parts) == 2:
            name, left, right = f"Target_{idx}", parts[0], parts[1]
        else:
            raise ValueError(f"Invalid target token: '{tok}'. Expected format: 'Name:LeftArm:RightArm'")
        targets.append({
            "name": name,
            "left_arm": left.upper(),
            "right_arm": right.upper(),
        })
    return targets


def main():
    parser = argparse.ArgumentParser(
        description="design_multiplex.py: Design orthogonal multiplex (triplex) padlock probe panels with shared RCA primer and non-cross-reactive detection probes."
    )
    parser.add_argument(
        "-j",
        "--job-name",
        required=True,
        help="Job name (outputs saved under results/<job_name>/)",
    )
    parser.add_argument(
        "-t",
        "--targets-csv",
        help="Path to CSV containing columns: Target_Name, Left_Arm, Right_Arm",
    )
    parser.add_argument(
        "--targets",
        help="Inline target specification: 'Name1:Left1:Right1,Name2:Left2:Right2,...'",
    )
    parser.add_argument(
        "--left-arm",
        default="",
        help="Left arm sequence (5'->3') for a new target",
    )
    parser.add_argument(
        "--right-arm",
        default="",
        help="Right arm sequence (5'->3') for a new target",
    )
    parser.add_argument(
        "--target-name",
        default="",
        help="Name for the new target (used with --left-arm and --right-arm, default: Target_New)",
    )

    # Fixed Reference Set options
    parser.add_argument(
        "--fixed-oligos",
        default="",
        help="Path to CSV file or inline string containing fixed reference oligos (PLP, RCA_Primer, PLP_Detection_Probe)",
    )
    parser.add_argument(
        "--fixed-plp",
        default="",
        help="Pre-existing fixed PLP sequence",
    )
    parser.add_argument(
        "--fixed-rca-primer",
        default="",
        help="Pre-existing fixed RCA primer sequence",
    )
    parser.add_argument(
        "--fixed-label",
        default="5'-Biotin",
        help="Reporter label for the fixed channel RCA primer (default: 5'-Biotin)",
    )

    parser.add_argument(
        "--universal-detection-probe",
        "--universal-probe",
        dest="universal_detection_probe",
        default="",
        help="Optional fixed universal detection probe sequence (5'->3', 5'-FAM labeled)",
    )
    parser.add_argument(
        "--universal-rca-primer",
        default="",
        help="Deprecated alias for fixed RCA primer (retained for backward compatibility)",
    )
    parser.add_argument(
        "--labels",
        default="5'-Biotin,5'-DIG,5'-TAMRA,5'-Cy5,5'-TexasRed",
        help="Comma-separated reporter labels for RCA primers (default: 5'-Biotin,5'-DIG,5'-TAMRA,5'-Cy5,5'-TexasRed)",
    )
    parser.add_argument(
        "--max-cross-dimer-dg",
        type=float,
        default=None,
        help="Maximum permissible cross-dimer binding dG in kcal/mol (default: from config or -2.0)",
    )
    parser.add_argument(
        "-c",
        "--config",
        default="config.json",
        help="Optional config.json to inherit thresholds from (default: config.json)",
    )
    parser.add_argument(
        "--restriction-site",
        default="",
        help="Optional restriction enzyme cleavage site (e.g. BstUI, EcoRV, HaeIII, or raw sequence like CGCG) inserted into PLP backbone",
    )

    args = parser.parse_args()

    # Determine job directory
    raw_job = str(args.job_name).strip()
    job_path = Path(raw_job)
    if not job_path.is_absolute():
        parts = job_path.parts
        if not parts or parts[0] != "results":
            job_path = Path("results") / job_path
    job_dir = str(job_path.resolve())

    # Load targets
    if args.targets_csv:
        targets = parse_targets_csv(args.targets_csv)
    elif args.targets:
        targets = parse_inline_targets(args.targets)
    elif args.left_arm and args.right_arm:
        t_name = args.target_name.strip() if args.target_name else "Target_New"
        targets = [{
            "name": t_name,
            "left_arm": args.left_arm.strip().upper(),
            "right_arm": args.right_arm.strip().upper(),
        }]
    else:
        print("[Error] Please specify targets via --targets-csv, --targets, or (--left-arm and --right-arm).", file=sys.stderr)
        sys.exit(1)

    # Determine fixed reference panel if provided
    fixed_panel = None
    if args.fixed_plp and args.fixed_rca_primer:
        fixed_panel = [{
            "plp": args.fixed_plp.strip().upper(),
            "rca_primer": args.fixed_rca_primer.strip().upper(),
            "detection_probe": args.universal_detection_probe.strip().upper(),
            "channel_label": args.fixed_label.strip(),
            "target_name": "Target_Fixed_1",
        }]
    elif args.fixed_oligos:
        fixed_panel = args.fixed_oligos

    labels = [lbl.strip() for lbl in args.labels.split(",") if lbl.strip()]

    # Load config if available
    cfg = {}
    if os.path.isfile(args.config):
        try:
            with open(args.config, "r") as f:
                cfg = json.load(f)
        except Exception:
            cfg = {}

    max_cross_dimer_dg = args.max_cross_dimer_dg
    if max_cross_dimer_dg is None:
        max_cross_dimer_dg = cfg.get("multiplex", {}).get("max_cross_dimer_dg", -2.0)

    try:
        design_multiplex_panel(
            targets=targets,
            output_dir=job_dir,
            universal_detection_probe=args.universal_detection_probe,
            universal_rca_primer=args.universal_rca_primer,
            channel_labels=labels,
            fixed_panel=fixed_panel,
            max_cross_dimer_dg=max_cross_dimer_dg,
            restriction_site=args.restriction_site,
            cfg=cfg,
        )
    except Exception as err:
        print(f"\n[Error] {err}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
