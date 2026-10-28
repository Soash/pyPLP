"""
pyPLP.crispr.designer: CRISPR Cas mode oligo and PAM designer.
Generates [A/C/G]-T-T-[A/C/G] (VTTV) PAM junction engineering, Cis probe, and crRNA guide sequences
for Padlock Probe (PLP) and Rolling Circle Amplification (RCA) detection.
"""

import re
from Bio.Seq import Seq

from pyPLP.thermodynamics.metrics import compute_metrics

DEFAULT_CRRNA_REPEAT = "UAAUUUCUACUAAGUGUAGAU"
DEFAULT_V = "A"

# Standard Cas12a collateral cleavage reporter probes
CRISPR_REPORTERS = [
    {
        "oligo_name": "CRISPR_FQ_Reporter_Probe",
        "sequence_5to3": "5'-FAM-TTATTATT-BHQ1-3'",
        "base_sequence": "TTATTATT",
        "modifications": "5'-FAM / 3'-BHQ1",
        "type": "DNA",
        "length_nt": 8,
        "gc_percent": 0.0,
        "tm_c": "N/A",
        "role": "CRISPR Cas12a Fluorescent Quenched Reporter (FQ Probe: 5'-FAM-TTATTATT-BHQ1-3')",
    },
    {
        "oligo_name": "CRISPR_FB_Reporter_Probe",
        "sequence_5to3": "5'-FAM-TTATTATT-Biotin-3'",
        "base_sequence": "TTATTATT",
        "modifications": "5'-FAM / 3'-Biotin",
        "type": "DNA",
        "length_nt": 8,
        "gc_percent": 0.0,
        "tm_c": "N/A",
        "role": "CRISPR Cas12a Lateral Flow Biosensor Reporter (FB Probe: 5'-FAM-TTATTATT-Biotin-3')",
    },
]


def determine_pam_and_rest_tail(right_arm: str, default_v: str = DEFAULT_V) -> dict:
    """
    Analyzes the 5' terminus and internal sequence of the Right Arm to determine
    the required Rest tail nucleotides to form an [A/C/G]-T-T-[A/C/G] (VTTV) PAM
    (V = A, C, G; not T) at the Rest-Right Arm junction or inside the Right Arm.
    """
    clean_ra = str(right_arm).strip().upper()
    if not clean_ra:
        raise ValueError("Right Arm sequence cannot be empty.")

    clean_v = default_v.upper() if default_v.upper() in ("A", "C", "G") else "A"

    # Case 1: Right Arm already starts with [ACG]TT[ACG]
    if re.match(r"^[ACG]TT[ACG]", clean_ra):
        return {
            "pam_mode": "ra_prefix",
            "rest_tail": "",
            "pam_seq": clean_ra[:4],
            "internal_pos": 0,
            "pam_motif": "[A/C/G]-T-T-[A/C/G]",
        }

    # Case 2: Right Arm starts with TT[ACG] -> Rest supplies V (1 nt)
    if re.match(r"^TT[ACG]", clean_ra):
        return {
            "pam_mode": "junction",
            "rest_tail": clean_v,
            "pam_seq": clean_v + clean_ra[:3],
            "internal_pos": None,
            "pam_motif": "[A/C/G]-T-T-[A/C/G]",
        }

    # Case 3: Right Arm starts with T[ACG] -> Rest supplies VT (2 nt)
    if re.match(r"^T[ACG]", clean_ra):
        return {
            "pam_mode": "junction",
            "rest_tail": clean_v + "T",
            "pam_seq": clean_v + "T" + clean_ra[:2],
            "internal_pos": None,
            "pam_motif": "[A/C/G]-T-T-[A/C/G]",
        }

    # Case 4: Right Arm starts with [ACG] -> Rest supplies VTT (3 nt)
    if re.match(r"^[ACG]", clean_ra):
        return {
            "pam_mode": "junction",
            "rest_tail": clean_v + "TT",
            "pam_seq": clean_v + "TT" + clean_ra[0],
            "internal_pos": None,
            "pam_motif": "[A/C/G]-T-T-[A/C/G]",
        }

    # Case 5: Right Arm starts with >= 3 Ts -> Check for internal [ACG]TT[ACG]
    match = re.search(r"[ACG]TT[ACG]", clean_ra)
    if match:
        return {
            "pam_mode": "internal_ra",
            "rest_tail": "",
            "pam_seq": match.group(0),
            "internal_pos": match.start(),
            "pam_motif": "[A/C/G]-T-T-[A/C/G]",
        }

    # Fallback: force full VTTV at junction from Rest
    return {
        "pam_mode": "junction_fallback",
        "rest_tail": clean_v + "TT" + clean_v,
        "pam_seq": clean_v + "TT" + clean_v,
        "internal_pos": None,
        "pam_motif": "[A/C/G]-T-T-[A/C/G]",
    }


def determine_shared_pam_and_rest_tail(targets: list[dict], default_v: str = DEFAULT_V) -> dict:
    """
    Analyzes Right Arm sequences across multiple targets to engineer a shared Rest tail
    that guarantees a valid Cas12a VTTV PAM ([A/C/G]-T-T-[A/C/G]) for every target probe.

    1. If all targets require the identical partial Rest tail (e.g. all need "A"), returns that tail.
    2. If targets have differing requirements or any target lacks a junction/internal PAM,
       returns a universal standalone VTTV PAM at Rest 3'-terminus ("CTTA").
    """
    if not targets:
        raise ValueError("Targets list cannot be empty.")

    clean_v = default_v.upper() if default_v.upper() in ("A", "C", "G") else "A"

    individual_info = []
    for t in targets:
        ra = t.get("right_arm", "")
        info = determine_pam_and_rest_tail(ra, default_v=clean_v)
        individual_info.append(info)

    tails = set(info["rest_tail"] for info in individual_info)
    modes = set(info["pam_mode"] for info in individual_info)

    # Check if all targets can be satisfied by the same exact partial tail
    if len(tails) == 1 and ("junction_fallback" not in modes):
        common_tail = individual_info[0]["rest_tail"]
        common_mode = individual_info[0]["pam_mode"]
        common_pam = individual_info[0]["pam_seq"]
        all_same_pam = (len(set(info["pam_seq"] for info in individual_info)) == 1)
        return {
            "pam_mode": f"shared_{common_mode}",
            "rest_tail": common_tail,
            "pam_seq": common_pam if all_same_pam else "Target-specific VTTV",
            "pam_location": f"Shared junction ({len(common_tail)} nt from Rest)" if common_tail else "Right Arm prefix/internal",
            "pam_motif": "[A/C/G]-T-T-[A/C/G]",
        }

    # Differing requirements or fallback needed across targets:
    # Use universal standalone VTTV PAM ("CTTA") at the Rest 3'-end.
    # This guarantees that every target Right Arm starts immediately 3' of the PAM.
    universal_pam = "CTTA"
    return {
        "pam_mode": "shared_universal_pam",
        "rest_tail": universal_pam,
        "pam_seq": universal_pam,
        "pam_location": f"Rest 3'-terminus (Universal Standalone VTTV PAM: {universal_pam})",
        "pam_motif": "[A/C/G]-T-T-[A/C/G]",
    }


def design_crispr_components(
    left_arm: str,
    right_arm: str,
    backbone: str = "",
    crrna_repeat: str = DEFAULT_CRRNA_REPEAT,
    spacer_len: int = 20,
    cis_probe_len: int = 30,
    generate_cis_probe: bool = False,
    default_v: str = DEFAULT_V,
) -> dict:
    """
    Assembles CRISPR Cas assay components for [A/C/G]-T-T-[A/C/G] (VTTV) PAM:
    1. Detects PAM position (junction vs internal).
    2. Constructs the Cis Probe (~25-31 nt DNA) spanning PAM across the circular ligation junction.
    3. Extracts the 20 nt protospacer and constructs the crRNA (fixed repeat + spacer).
    """
    clean_la = str(left_arm).strip().upper()
    clean_ra = str(right_arm).strip().upper()
    clean_bb = str(backbone).strip().upper() if backbone else ""

    if clean_bb:
        bb_tail = clean_bb[-4:] if len(clean_bb) >= 4 else clean_bb
        full_junction = bb_tail + clean_ra + clean_la + clean_bb
        ra_start = len(bb_tail)
        ra_end = ra_start + len(clean_ra)
        target_span_end = ra_end + len(clean_la)

        # Search for [ACG]TT[ACG] starting strictly within bb_tail or right_arm (upstream of ligation nick)
        match = None
        for m in re.finditer(r"[ACG]TT[ACG]", full_junction):
            if m.start() < ra_end:
                match = m
                break

        if match:
            start = match.start()
            pam_seq = match.group(0)
            if start < ra_start:
                pam_mode = "junction"
                bases_from_bb = ra_start - start
                pam_location_desc = f"Junction (Rest: +{bases_from_bb} nt, RA: {4 - bases_from_bb} nt)"
                rest_tail = clean_bb[len(clean_bb) - bases_from_bb :]
            else:
                pam_mode = "internal_ra"
                pos = start - ra_start
                pam_location_desc = f"Internal Right Arm (pos {pos + 1})"
                rest_tail = ""

            junction_seq = full_junction[start:]
            protospacer_dna = junction_seq[4 : 4 + spacer_len]
            rem_arm_len = target_span_end - start
            if 25 <= rem_arm_len <= max(31, cis_probe_len):
                cis_probe_seq = junction_seq[:rem_arm_len]
            else:
                cis_probe_seq = junction_seq[:cis_probe_len]
            valid_pam = True
        else:
            # No valid PAM upstream of ligation nick
            pam_mode = "none"
            rest_tail = ""
            pam_seq = "None"
            pam_location_desc = "None (No valid VTTV PAM in Rest tail or Right Arm)"
            protospacer_dna = ""
            cis_probe_seq = ""
            valid_pam = False
    else:
        pam_info = determine_pam_and_rest_tail(clean_ra, default_v=default_v)
        pam_mode = pam_info["pam_mode"]
        rest_tail = pam_info["rest_tail"]
        pam_seq = pam_info["pam_seq"]
        valid_pam = (pam_mode != "none")

        if pam_mode in ("junction", "junction_fallback"):
            junction_seq = rest_tail + clean_ra + clean_la + clean_bb
            protospacer_dna = junction_seq[4 : 4 + spacer_len]
            pam_location_desc = f"Junction (Rest: +{len(rest_tail)} nt, RA: {4 - len(rest_tail)} nt)"
            rem_arm_len = len(rest_tail) + len(clean_ra) + len(clean_la)
            cis_probe_seq = junction_seq[:rem_arm_len] if 25 <= rem_arm_len <= max(31, cis_probe_len) else junction_seq[:cis_probe_len]
        elif pam_mode == "internal_ra":
            pos = pam_info["internal_pos"]
            junction_seq = clean_ra[pos:] + clean_la + clean_bb
            protospacer_dna = junction_seq[4 : 4 + spacer_len]
            pam_location_desc = f"Internal Right Arm (pos {pos + 1})"
            rem_arm_len = (len(clean_ra) - pos) + len(clean_la)
            cis_probe_seq = junction_seq[:rem_arm_len] if 25 <= rem_arm_len <= max(31, cis_probe_len) else junction_seq[:cis_probe_len]
        else:  # ra_prefix
            junction_seq = clean_ra + clean_la + clean_bb
            protospacer_dna = junction_seq[4 : 4 + spacer_len]
            pam_location_desc = "Right Arm 5'-terminus"
            rem_arm_len = len(clean_ra) + len(clean_la)
            cis_probe_seq = junction_seq[:rem_arm_len] if 25 <= rem_arm_len <= max(31, cis_probe_len) else junction_seq[:cis_probe_len]

    # Clean crRNA repeat (ensure RNA format)
    clean_repeat = crrna_repeat.strip().upper().replace("T", "U")

    if valid_pam and protospacer_dna:
        # Convert protospacer DNA to RNA spacer
        spacer_rna = protospacer_dna.replace("T", "U")
        full_crrna = f"{clean_repeat}{spacer_rna}"
        crrna_display = f"{clean_repeat} {spacer_rna}"
        spacer_gc, spacer_tm = compute_metrics(Seq(protospacer_dna))
    else:
        spacer_rna = ""
        full_crrna = ""
        crrna_display = "None"
        spacer_gc, spacer_tm = 0.0, 0.0

    # Calculate biophysical metrics
    if generate_cis_probe and valid_pam and cis_probe_seq:
        cis_gc, cis_tm = compute_metrics(Seq(cis_probe_seq))
        cis_len = len(cis_probe_seq)
    else:
        cis_probe_seq = ""
        cis_len = 0
        cis_gc = 0.0
        cis_tm = 0.0

    return {
        "valid_pam": valid_pam,
        "pam_seq": pam_seq,
        "pam_mode": pam_mode,
        "pam_location": pam_location_desc,
        "rest_tail": rest_tail,
        "cis_probe_seq": cis_probe_seq,
        "cis_probe_len": cis_len,
        "cis_probe_gc": cis_gc,
        "cis_probe_tm": cis_tm,
        "crrna_repeat": clean_repeat,
        "spacer_seq_rna": spacer_rna,
        "spacer_seq_dna": protospacer_dna,
        "spacer_len": len(spacer_rna),
        "spacer_gc": spacer_gc,
        "spacer_tm": spacer_tm,
        "crrna_full_seq": full_crrna,
        "crrna_full_display": crrna_display,
        "crrna_full_len": len(full_crrna),
    }

