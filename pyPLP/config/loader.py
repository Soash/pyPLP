import json
import os
import sys


def load_config(config_path: str = "config.json") -> dict:
    """Loads and validates pipeline configuration.

    Supports both single arm sizes and range intervals (e.g., 12 or [12, 20]).
    """
    if not os.path.exists(config_path):
        print(f"Error: Config file '{config_path}' not found.", file=sys.stderr)
        sys.exit(1)

    with open(config_path, "r") as f:
        cfg = json.load(f)

    fasta_file = cfg.get("input", "input.fasta")

    # Parse arm size: can be single int/str (e.g., 18, "18") or a range list [12, 20]
    raw_arm_size = cfg.get("armSize", 12)
    if isinstance(raw_arm_size, list):
        if len(raw_arm_size) == 2:
            start_sz, end_sz = int(raw_arm_size[0]), int(raw_arm_size[1])
            arm_sizes = list(range(start_sz, end_sz + 1))
        else:
            arm_sizes = [int(x) for x in raw_arm_size]
    elif isinstance(raw_arm_size, (int, str)) and str(raw_arm_size).isdigit():
        arm_sizes = [int(raw_arm_size)]
    else:
        print(
            f"Warning: Unrecognized armSize format '{raw_arm_size}'. Defaulting to [12]."
        )
        arm_sizes = [12]

    min_arm_tm_val = cfg.get("minArmTm")
    min_arm_tm = (
        float(min_arm_tm_val)
        if min_arm_tm_val not in (None, "")
        else None
    )

    hp_filter_val = cfg.get("Hairpin_Found")
    hairpin_filter = (
        str(hp_filter_val).strip().capitalize()
        if hp_filter_val not in (None, "")
        else None
    )

    max_homo_val = cfg.get("maxHomopolymer")
    max_homo = (
        int(max_homo_val)
        if max_homo_val not in (None, "")
        else None
    )

    blastdb_dir = cfg.get("blastdbDir", "blastdb")
    threads = int(cfg.get("threads", os.cpu_count() or 32))
    blast_min_ident = float(cfg.get("blastMinIdentity", 70.0))
    max_allowed_hits = cfg.get("maxAllowedHits")

    return {
        "fasta_file": fasta_file,
        "arm_sizes": arm_sizes,
        "min_arm_tm": min_arm_tm,
        "hairpin_filter": hairpin_filter,
        "max_homo": max_homo,
        "blastdb_dir": blastdb_dir,
        "threads": threads,
        "blast_min_ident": blast_min_ident,
        "max_allowed_hits": max_allowed_hits,
        "raw_config": cfg,
    }
