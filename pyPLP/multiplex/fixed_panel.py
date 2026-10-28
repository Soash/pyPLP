"""pyPLP.multiplex.fixed_panel: Parsing and biophysical validation of pre-existing fixed PLP panels.
"""

import csv
import io
import os
from Bio.Seq import Seq

from pyPLP.multiplex.cross_dimer import calc_heterodimer_dg, calc_homodimer_dg
from pyPLP.plp_generation.structure_eval import eval_plp_structure
from pyPLP.thermodynamics.hairpin import check_hairpin
from pyPLP.thermodynamics.homopolymer import get_max_homopolymer
from pyPLP.thermodynamics.metrics import compute_metrics
from pyPLP.thermodynamics.tm import calculate_tm


def parse_fixed_oligos(source: str | list | dict) -> dict:
    """Parses a fixed reference set of oligonucleotides (PLP, RCA Primer, Detection Probe).

    Supports:
    1. File path to CSV.
    2. Raw multiline CSV string (e.g., 'PLP,ATACC...\\nRCA_Primer,CGGTT...').
    3. Dict with keys 'plp', 'rca_primer', 'detection_probe' (or variations).
    4. List of dicts (for multiple fixed channels).

    Returns a normalized dictionary:
    {
        "universal_detection_probe": str,
        "channels": [
            {
                "target_name": str,
                "plp_seq": str,
                "rca_primer_seq": str,
                "detection_probe_seq": str,
                "channel_label": str,
            }
        ]
    }
    """
    if not source:
        return {"universal_detection_probe": "", "channels": []}

    channels = []
    univ_probe = ""

    # Case A: Dict
    if isinstance(source, dict):
        if "channels" in source or "fixed_panel" in source:
            sub = source.get("channels") or source.get("fixed_panel") or []
            return parse_fixed_oligos(sub)
        p_seq = source.get("plp") or source.get("plp_seq") or source.get("PLP") or ""
        rca_seq = source.get("rca_primer") or source.get("rca_primer_seq") or source.get("RCA_Primer") or ""
        probe_seq = source.get("detection_probe") or source.get("detection_probe_seq") or source.get("PLP_Detection_Probe") or source.get("universal_probe") or ""
        label = source.get("channel_label") or source.get("label") or "5'-Biotin"
        name = source.get("target_name") or source.get("name") or "Fixed_Target_1"
        if p_seq and rca_seq:
            channels.append({
                "target_name": str(name).strip(),
                "plp_seq": str(p_seq).strip().upper(),
                "rca_primer_seq": str(rca_seq).strip().upper(),
                "detection_probe_seq": str(probe_seq).strip().upper(),
                "channel_label": str(label).strip(),
            })
            if probe_seq:
                univ_probe = str(probe_seq).strip().upper()
        return {"universal_detection_probe": univ_probe, "channels": channels}

    # Case B: List of dicts
    if isinstance(source, list):
        for idx, item in enumerate(source, start=1):
            if isinstance(item, dict):
                p_seq = item.get("plp") or item.get("plp_seq") or item.get("PLP") or ""
                rca_seq = item.get("rca_primer") or item.get("rca_primer_seq") or item.get("RCA_Primer") or ""
                probe_seq = item.get("detection_probe") or item.get("detection_probe_seq") or item.get("PLP_Detection_Probe") or item.get("universal_probe") or ""
                label = item.get("channel_label") or item.get("label") or f"5'-Biotin"
                name = item.get("target_name") or item.get("name") or f"Fixed_Target_{idx}"
                if p_seq and rca_seq:
                    channels.append({
                        "target_name": str(name).strip(),
                        "plp_seq": str(p_seq).strip().upper(),
                        "rca_primer_seq": str(rca_seq).strip().upper(),
                        "detection_probe_seq": str(probe_seq).strip().upper(),
                        "channel_label": str(label).strip(),
                    })
                    if probe_seq and not univ_probe:
                        univ_probe = str(probe_seq).strip().upper()
        return {"universal_detection_probe": univ_probe, "channels": channels}

    # Case C: String (file path or multiline CSV text)
    if isinstance(source, str):
        content = ""
        if os.path.isfile(source):
            with open(source, "r", encoding="utf-8") as f:
                content = f.read()
        else:
            content = source

        lines = [line.strip() for line in content.splitlines() if line.strip()]
        if not lines:
            return {"universal_detection_probe": "", "channels": []}

        # Check if tabular with header
        first_line = lines[0].lower()
        if ("target_name" in first_line or "plp_sequence" in first_line or "rca_primer" in first_line) and "," in lines[0] and len(lines[0].split(",")) > 2:
            # Tabular CSV
            reader = csv.DictReader(io.StringIO(content))
            for idx, row in enumerate(reader, start=1):
                p_seq = row.get("PLP_Sequence") or row.get("PLP") or row.get("plp_seq") or ""
                rca_seq = row.get("RCA_Primer_Sequence") or row.get("RCA_Primer") or row.get("rca_primer") or ""
                probe_seq = row.get("Detection_Probe_Sequence") or row.get("Detection_Probe") or row.get("PLP_Detection_Probe") or row.get("universal_probe") or ""
                name = row.get("Target_Name") or row.get("Name") or f"Fixed_Target_{idx}"
                label = row.get("Channel_Label") or row.get("Label") or "5'-Biotin"
                if p_seq and rca_seq:
                    channels.append({
                        "target_name": str(name).strip(),
                        "plp_seq": str(p_seq).strip().upper(),
                        "rca_primer_seq": str(rca_seq).strip().upper(),
                        "detection_probe_seq": str(probe_seq).strip().upper(),
                        "channel_label": str(label).strip(),
                    })
                    if probe_seq and not univ_probe:
                        univ_probe = str(probe_seq).strip().upper()
            return {"universal_detection_probe": univ_probe, "channels": channels}

        # Simple Key-Value CSV lines (e.g., PLP,seq / RCA_Primer,seq / PLP_Detection_Probe,seq)
        current_target = "Fixed_Target_1"
        curr_plp = ""
        curr_rca = ""
        curr_probe = ""
        curr_label = "5'-Biotin"

        for line in lines:
            parts = [p.strip() for p in line.split(",") if p.strip()]
            if len(parts) < 2:
                continue
            key = parts[0].strip().upper()
            val = parts[1].strip().upper()

            if key in ("PLP", "PLP_SEQUENCE", "PADLOCK", "PADLOCK_PROBE"):
                curr_plp = val
            elif key in ("RCA_PRIMER", "RCA", "PRIMER", "RCA_PRIMER_SEQUENCE"):
                curr_rca = val
                if len(parts) >= 3 and parts[2].strip():
                    curr_label = parts[2].strip()
            elif key in ("PLP_DETECTION_PROBE", "DETECTION_PROBE", "PROBE", "UNIVERSAL_PROBE"):
                curr_probe = val
                univ_probe = val
            elif key in ("TARGET", "TARGET_NAME", "NAME"):
                current_target = parts[1].strip()

        if curr_plp and curr_rca:
            channels.append({
                "target_name": current_target,
                "plp_seq": curr_plp,
                "rca_primer_seq": curr_rca,
                "detection_probe_seq": curr_probe or univ_probe,
                "channel_label": curr_label,
            })

    return {"universal_detection_probe": univ_probe, "channels": channels}


def evaluate_fixed_channel(ch_data: dict, channel_idx: int = 1) -> dict:
    """Evaluates the thermodynamic secondary structure and biophysical properties of a fixed PLP channel."""
    plp_seq = ch_data["plp_seq"]
    rca_seq = ch_data["rca_primer_seq"]
    probe_seq = ch_data.get("detection_probe_seq", "")
    target_name = ch_data.get("target_name", f"Fixed_Target_{channel_idx}")
    channel_label = ch_data.get("channel_label", "5'-Biotin")

    # Evaluate PLP
    full_gc, full_tm = compute_metrics(Seq(plp_seq))
    dot_bracket, mfe_dg = eval_plp_structure(plp_seq)
    plp_id = f"PLP_{channel_idx}_{target_name}_Fixed"

    # Evaluate RCA Primer
    rca_tm = calculate_tm(rca_seq)
    rca_gc, _ = compute_metrics(Seq(rca_seq))
    rca_homo_len, rca_homo_seq = get_max_homopolymer(rca_seq)
    _, rca_hp_tm, _ = check_hairpin(rca_seq)
    rca_self_dg = calc_homodimer_dg(rca_seq)
    rca_binding_dg = calc_heterodimer_dg(plp_seq, rca_seq)

    # Evaluate Detection probe
    probe_tm = calculate_tm(probe_seq) if probe_seq else 0.0
    probe_gc = compute_metrics(Seq(probe_seq))[0] if probe_seq else 0.0

    return {
        "plp_id": plp_id,
        "target_name": target_name,
        "channel_idx": channel_idx,
        "channel_label": channel_label,
        "full_plp": plp_seq,
        "total_length": len(plp_seq),
        "full_gc": full_gc,
        "full_tm": full_tm,
        "dot_bracket": dot_bracket,
        "mfe_dg": mfe_dg,
        "is_fixed": True,
        "rca_primer_seq": rca_seq,
        "rca_primer_tm": rca_tm,
        "rca_primer_gc": rca_gc,
        "rca_primer_homo": f"{rca_homo_len} ({rca_homo_seq})",
        "rca_primer_hp_tm": rca_hp_tm,
        "rca_primer_self_dg": rca_self_dg,
        "rca_primer_binding_dg": rca_binding_dg,
        "rca_primer_mod": channel_label,
        "detection_probe_seq": probe_seq,
        "detection_probe_tm": probe_tm,
        "detection_probe_gc": probe_gc,
        "detection_probe_mod": "5'-FAM",
    }
