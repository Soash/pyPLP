"""pyPLP.multiplex.probe_generator: Automated design of mutually orthogonal detection probes.
"""

import random
from Bio.Seq import Seq

from pyPLP.multiplex.cross_dimer import calc_heterodimer_dg, calc_homodimer_dg
from pyPLP.thermodynamics.hairpin import check_hairpin
from pyPLP.thermodynamics.homopolymer import get_max_homopolymer
from pyPLP.thermodynamics.metrics import compute_metrics
from pyPLP.thermodynamics.tm import calculate_tm


def generate_candidate_probe(
    length: int = 20,
    target_tm: float = 50.0,
    tm_tol: float = 3.0,
    max_homo: int = 4,
    max_hairpin_tm: float = 30.0,
    max_attempts: int = 500,
) -> dict | None:
    """Generates a single candidate detection probe meeting biophysical criteria."""
    for _ in range(max_attempts):
        seq = "".join(random.choices(["A", "C", "G", "T"], k=length))

        homo_len, homo_seq = get_max_homopolymer(seq)
        if homo_len >= max_homo:
            continue

        hp_found, hp_tm, _ = check_hairpin(seq)
        if hp_found == "Yes" and hp_tm > max_hairpin_tm:
            continue

        tm = calculate_tm(seq)
        if abs(tm - target_tm) <= tm_tol:
            self_dg = calc_homodimer_dg(seq)
            gc, _ = compute_metrics(Seq(seq))
            return {
                "seq": seq,
                "tm": tm,
                "gc": gc,
                "homo": f"{homo_len} ({homo_seq})",
                "hp_tm": hp_tm,
                "self_dg": self_dg,
            }

    return None


def generate_universal_detection_probe(
    probe_length: int = 15,
    target_tm: float = 50.0,
    tm_tol: float = 3.5,
    max_cross_dg: float = -2.0,
    max_homo: int = 4,
    max_hairpin_tm: float = 28.0,
    predefined_seq: str = "",
    target_arms_list: list[tuple[str, str]] | None = None,
    max_attempts: int = 3000,
) -> dict | None:
    """Generates a universal detection probe (5'-FAM labeled) compatible with all target arms."""
    if predefined_seq:
        seq = str(predefined_seq).strip().upper()
        homo_len, homo_seq = get_max_homopolymer(seq)
        _, hp_tm, _ = check_hairpin(seq)
        tm = calculate_tm(seq)
        gc, _ = compute_metrics(Seq(seq))
        self_dg = calc_homodimer_dg(seq)
        return {
            "seq": seq,
            "tm": tm,
            "gc": gc,
            "homo": f"{homo_len} ({homo_seq})",
            "hp_tm": hp_tm,
            "self_dg": self_dg,
            "label": "5'-FAM",
        }

    for _ in range(max_attempts):
        cand = generate_candidate_probe(
            length=probe_length,
            target_tm=target_tm,
            tm_tol=tm_tol,
            max_homo=max_homo,
            max_hairpin_tm=max_hairpin_tm,
        )
        if not cand:
            continue
        if cand["self_dg"] < max_cross_dg:
            continue

        clash = False
        if target_arms_list:
            for l_arm, r_arm in target_arms_list:
                if (
                    calc_heterodimer_dg(cand["seq"], l_arm) < max_cross_dg
                    or calc_heterodimer_dg(cand["seq"], r_arm) < max_cross_dg
                ):
                    clash = True
                    break
        if clash:
            continue

        cand["label"] = "5'-FAM"
        return cand

    return None


def generate_orthogonal_rca_primers(
    channel_count: int = 3,
    universal_probe_seq: str = "",
    primer_length: int = 15,
    target_tm: float = 50.0,
    tm_tol: float = 3.5,
    max_cross_dg: float = -2.0,
    max_homo: int = 4,
    max_hairpin_tm: float = 28.0,
    predefined_primers: list[str] | None = None,
    channel_labels: list[str] | None = None,
    fixed_plps: list[str] | None = None,
    fixed_rca_primers: list[str] | None = None,
    target_arms_list: list[tuple[str, str]] | None = None,
    max_search_cycles: int = 3000,
) -> list[dict] | None:
    """Generates a mutually orthogonal set of K channel-specific RCA primers compatible with
    a shared universal detection probe, free of inter-primer, primer-arm, and fixed-panel cross-reactivity.
    """
    default_labels = ["5'-Biotin", "5'-DIG", "5'-TAMRA", "5'-Cy5", "5'-TexasRed"]
    labels = channel_labels or default_labels[:channel_count]
    if len(labels) < channel_count:
        labels.extend([f"Channel_{i+1}" for i in range(len(labels), channel_count)])

    probe_seq = str(universal_probe_seq).strip().upper() if universal_probe_seq else ""

    for _ in range(max_search_cycles):
        selected_primers = []
        success = True

        for ch_idx in range(channel_count):
            label = labels[ch_idx]

            # Use predefined primer if provided for this channel
            if predefined_primers and ch_idx < len(predefined_primers) and predefined_primers[ch_idx]:
                p_seq = str(predefined_primers[ch_idx]).strip().upper()
                p_tm = calculate_tm(p_seq)
                p_gc, _ = compute_metrics(Seq(p_seq))
                p_homo_len, p_homo_seq = get_max_homopolymer(p_seq)
                _, p_hp_tm, _ = check_hairpin(p_seq)
                cand = {
                    "seq": p_seq,
                    "tm": p_tm,
                    "gc": p_gc,
                    "homo": f"{p_homo_len} ({p_homo_seq})",
                    "hp_tm": p_hp_tm,
                    "self_dg": calc_homodimer_dg(p_seq),
                }
            else:
                cand = generate_candidate_probe(
                    length=primer_length,
                    target_tm=target_tm,
                    tm_tol=tm_tol,
                    max_homo=max_homo,
                    max_hairpin_tm=max_hairpin_tm,
                )

            if not cand:
                success = False
                break

            # Check self-homodimer
            if cand["self_dg"] < max_cross_dg:
                success = False
                break

            cand_rc = str(Seq(cand["seq"]).reverse_complement())

            # Check orthogonality against universal detection probe (both sense and rc for backbone)
            if probe_seq:
                probe_cross_dg = calc_heterodimer_dg(cand["seq"], probe_seq)
                probe_rc_dg = calc_heterodimer_dg(cand_rc, probe_seq)
                if probe_cross_dg < max_cross_dg or probe_rc_dg < max_cross_dg:
                    success = False
                    break

            # Check orthogonality against fixed PLPs and fixed RCA primers (only for newly generated primers)
            is_predefined = bool(predefined_primers and ch_idx < len(predefined_primers) and predefined_primers[ch_idx])
            if not is_predefined:
                if fixed_rca_primers:
                    clash_fixed_rca = False
                    for f_rca in fixed_rca_primers:
                        if f_rca and (
                            calc_heterodimer_dg(cand["seq"], f_rca) < max_cross_dg
                            or calc_heterodimer_dg(cand_rc, f_rca) < max_cross_dg
                        ):
                            clash_fixed_rca = True
                            break
                    if clash_fixed_rca:
                        success = False
                        break

                if fixed_plps:
                    clash_fixed_plp = False
                    for f_plp in fixed_plps:
                        if f_plp and (
                            calc_heterodimer_dg(cand["seq"], f_plp) < max_cross_dg
                            or calc_heterodimer_dg(cand_rc, f_plp) < max_cross_dg
                        ):
                            clash_fixed_plp = True
                            break
                    if clash_fixed_plp:
                        success = False
                        break

            # Check orthogonality against all target arms to prevent primer-arm cross-hybridization
            if target_arms_list:
                clash_arm = False
                for l_arm, r_arm in target_arms_list:
                    if (
                        calc_heterodimer_dg(cand["seq"], l_arm) < max_cross_dg
                        or calc_heterodimer_dg(cand["seq"], r_arm) < max_cross_dg
                        or calc_heterodimer_dg(cand_rc, l_arm) < max_cross_dg
                        or calc_heterodimer_dg(cand_rc, r_arm) < max_cross_dg
                    ):
                        clash_arm = True
                        break
                if clash_arm:
                    success = False
                    break

            # Check orthogonality against already selected RCA primers
            clash = False
            for prev in selected_primers:
                pair_dg = calc_heterodimer_dg(cand["seq"], prev["seq"])
                if pair_dg < max_cross_dg:
                    clash = True
                    break

            if clash:
                success = False
                break

            cand["channel_idx"] = ch_idx + 1
            cand["label"] = label
            selected_primers.append(cand)

        if success and len(selected_primers) == channel_count:
            return selected_primers

    return None


def generate_orthogonal_probe_set(
    channel_count: int = 3,
    rca_primer_seq: str = "",
    probe_length: int = 20,
    target_tm: float = 50.0,
    tm_tol: float = 3.0,
    max_cross_dg: float = -2.0,
    max_homo: int = 4,
    max_hairpin_tm: float = 30.0,
    predefined_probes: list[str] | None = None,
    channel_labels: list[str] | None = None,
    target_arms_list: list[tuple[str, str]] | None = None,
    max_search_cycles: int = 3000,
) -> list[dict] | None:
    """Generates a mutually orthogonal set of K detection probes compatible with
    a shared RCA primer and free of inter-probe and probe-arm cross-reactivity.
    """
    labels = channel_labels or [f"Channel_{i+1}" for i in range(channel_count)]
    if len(labels) < channel_count:
        labels.extend([f"Channel_{i+1}" for i in range(len(labels), channel_count)])

    rca_seq = str(rca_primer_seq).strip().upper() if rca_primer_seq else ""

    for _ in range(max_search_cycles):
        selected_probes = []
        success = True

        for ch_idx in range(channel_count):
            label = labels[ch_idx]

            # Use predefined probe if provided for this channel
            if predefined_probes and ch_idx < len(predefined_probes) and predefined_probes[ch_idx]:
                p_seq = str(predefined_probes[ch_idx]).strip().upper()
                p_tm = calculate_tm(p_seq)
                p_gc, _ = compute_metrics(Seq(p_seq))
                p_homo_len, p_homo_seq = get_max_homopolymer(p_seq)
                _, p_hp_tm, _ = check_hairpin(p_seq)
                cand = {
                    "seq": p_seq,
                    "tm": p_tm,
                    "gc": p_gc,
                    "homo": f"{p_homo_len} ({p_homo_seq})",
                    "hp_tm": p_hp_tm,
                    "self_dg": calc_homodimer_dg(p_seq),
                }
            else:
                cand = generate_candidate_probe(
                    length=probe_length,
                    target_tm=target_tm,
                    tm_tol=tm_tol,
                    max_homo=max_homo,
                    max_hairpin_tm=max_hairpin_tm,
                )

            if not cand:
                success = False
                break

            # Check self-homodimer
            if cand["self_dg"] < max_cross_dg:
                success = False
                break

            # Check orthogonality against universal RCA primer
            if rca_seq:
                rca_cross_dg = calc_heterodimer_dg(cand["seq"], rca_seq)
                if rca_cross_dg < max_cross_dg:
                    success = False
                    break

            # Check orthogonality against channel's own capture arms to prevent internal stems
            if target_arms_list and ch_idx < len(target_arms_list):
                l_arm, r_arm = target_arms_list[ch_idx]
                if calc_heterodimer_dg(cand["seq"], l_arm) < max_cross_dg:
                    success = False
                    break
                if calc_heterodimer_dg(cand["seq"], r_arm) < max_cross_dg:
                    success = False
                    break

            # Check orthogonality against already selected probes
            clash = False
            for prev in selected_probes:
                pair_dg = calc_heterodimer_dg(cand["seq"], prev["seq"])
                if pair_dg < max_cross_dg:
                    clash = True
                    break

            if clash:
                success = False
                break

            cand["channel_idx"] = ch_idx + 1
            cand["label"] = label
            selected_probes.append(cand)

        if success and len(selected_probes) == channel_count:
            return selected_probes

    return None
