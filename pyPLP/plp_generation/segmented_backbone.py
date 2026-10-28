import random
from Bio.Seq import Seq

from pyPLP.thermodynamics.hairpin import check_hairpin
from pyPLP.thermodynamics.homopolymer import get_max_homopolymer
from pyPLP.thermodynamics.metrics import compute_metrics
from pyPLP.thermodynamics.tm import calculate_tm


def generate_random_oligo_matching_tm(
    length: int,
    min_tm: float = 40.0,
    max_homo: int = 4,
    max_hairpin_tm: float = 35.0,
    max_attempts: int = 300,
) -> str:
    """Generates a random nucleotide sequence of specified length satisfying Tm,
    homopolymer, and hairpin constraints.
    """
    best_seq = None
    best_tm = -999.0

    for _ in range(max_attempts):
        candidate = "".join(random.choices(["A", "C", "G", "T"], k=length))

        homo_len, _ = get_max_homopolymer(candidate)
        if homo_len >= max_homo:
            continue

        hp_found, hp_tm, _ = check_hairpin(candidate)
        if hp_found == "Yes" and hp_tm > max_hairpin_tm:
            continue

        tm = calculate_tm(candidate)
        if tm >= min_tm:
            return candidate

        if tm > best_tm:
            best_tm = tm
            best_seq = candidate

    return best_seq or candidate


COMMON_RESTRICTION_ENZYMES = {
    "BSTUI": "CGCG",
    "ECORV": "GATATC",
    "HAEIII": "GGCC",
    "ALUI": "AGCT",
    "RSAI": "GTAC",
    "DPNI": "GATC",
    "BAMHI": "GGATCC",
    "PSTI": "CTGCAG",
    "HINDIII": "AAGCTT",
    "ECORI": "GAATTC",
    "NOTI": "GCGGCCGC",
    "XHOI": "CTCGAG",
    "SALGI": "GTCGAC",
    "SMAI": "CCCGGG",
}


def resolve_restriction_site(enzyme_or_seq: str) -> str:
    """Translates a restriction enzyme name (e.g. BstUI, EcoRV) or validates
    a raw DNA recognition sequence (e.g. CGCG, GATATC).
    """
    raw = str(enzyme_or_seq).strip().upper()
    if not raw:
        return ""
    if raw in COMMON_RESTRICTION_ENZYMES:
        return COMMON_RESTRICTION_ENZYMES[raw]
    if all(c in "ACGT" for c in raw):
        return raw
    raise ValueError(
        f"Invalid restriction enzyme or recognition sequence: '{enzyme_or_seq}'. "
        f"Supported enzymes: {', '.join(sorted(COMMON_RESTRICTION_ENZYMES.keys()))} or custom ACGT sequence."
    )


def build_initial_segmented_backbone(
    spacer1_len: int = 4,
    rca_primer_binding_len: int = 20,
    spacer2_len: int = 4,
    probe_binding_len: int = 20,
    rest_len: int = 4,
    restriction_site: str = "",
    universal_rca_primer_seq: str = "",
    universal_probe_seq: str = "",
    min_tm: float = 40.0,
    max_homo: int = 4,
    max_hairpin_tm: float = 35.0,
    spacer3_len: int = 0,
    re_mode: str = "standalone",
    re_probe_offset: int = 2,
    crispr_mode: bool = False,
    rest_tail: str = "",
) -> tuple[str, set[int]]:
    """Builds an initial segmented backbone composed of:
    Spacer1 + RCA_Primer_Binding_Region + Spacer2 + Probe_Binding_Region + Spacer3 + Restriction_Site + Rest
    Returns (backbone_str, protected_bb_indices).
    """
    # Spacer 1
    spacer1 = "".join(random.choices(["A", "C", "G", "T"], k=spacer1_len))
    while get_max_homopolymer(spacer1)[0] >= max_homo:
        spacer1 = "".join(random.choices(["A", "C", "G", "T"], k=spacer1_len))

    # RCA Primer Binding Region
    protected_indices = set()
    rca_universal = universal_rca_primer_seq.strip().upper() if universal_rca_primer_seq else ""
    if rca_universal:
        rca_binding = str(Seq(rca_universal).reverse_complement())
        start_idx = len(spacer1)
        for i in range(start_idx, start_idx + len(rca_binding)):
            protected_indices.add(i)
    else:
        primer_seq = generate_random_oligo_matching_tm(
            rca_primer_binding_len, min_tm=min_tm, max_homo=max_homo, max_hairpin_tm=max_hairpin_tm
        )
        rca_binding = str(Seq(primer_seq).reverse_complement())

    # Spacer 2
    spacer2 = "".join(random.choices(["A", "C", "G", "T"], k=spacer2_len))
    while get_max_homopolymer(spacer2)[0] >= max_homo:
        spacer2 = "".join(random.choices(["A", "C", "G", "T"], k=spacer2_len))

    re_seq = resolve_restriction_site(restriction_site)

    # Probe Binding Region
    probe_universal = universal_probe_seq.strip().upper() if universal_probe_seq else ""
    if probe_universal:
        probe_binding = probe_universal
        start_idx = len(spacer1) + len(rca_binding) + len(spacer2)
        for i in range(start_idx, start_idx + len(probe_binding)):
            protected_indices.add(i)
    else:
        if re_mode == "embedded_in_probe" and re_seq:
            # Generate probe binding with embedded restriction site
            right_pad_len = probe_binding_len - re_probe_offset - len(re_seq)
            if right_pad_len < 0:
                raise ValueError("probe_binding_len too small to embed restriction site at given offset.")
            while True:
                left_pad = "".join(random.choices(["A", "C", "G", "T"], k=re_probe_offset))
                right_pad = "".join(random.choices(["A", "C", "G", "T"], k=right_pad_len))
                cand = left_pad + re_seq + right_pad
                cand_probe = str(Seq(cand).reverse_complement())
                _, tm = compute_metrics(Seq(cand_probe))
                if tm >= min_tm and get_max_homopolymer(cand)[0] < max_homo:
                    _, hp_tm, _ = check_hairpin(cand_probe)
                    if hp_tm <= max_hairpin_tm:
                        probe_binding = cand
                        break
            # Protect the RE site within the probe binding region
            re_start = len(spacer1) + len(rca_binding) + len(spacer2) + re_probe_offset
            for i in range(re_start, re_start + len(re_seq)):
                protected_indices.add(i)
        else:
            probe_binding = generate_random_oligo_matching_tm(
                probe_binding_len, min_tm=min_tm, max_homo=max_homo, max_hairpin_tm=max_hairpin_tm
            )

    # Spacer 3
    spacer3 = "".join(random.choices(["A", "C", "G", "T"], k=spacer3_len)) if spacer3_len > 0 else ""
    while get_max_homopolymer(spacer3)[0] >= max_homo if spacer3 else False:
        spacer3 = "".join(random.choices(["A", "C", "G", "T"], k=spacer3_len))

    # Restriction Site (optional, standalone)
    standalone_re_seq = re_seq if re_mode == "standalone" else ""
    if standalone_re_seq:
        re_start = len(spacer1) + len(rca_binding) + len(spacer2) + len(probe_binding) + len(spacer3)
        for i in range(re_start, re_start + len(standalone_re_seq)):
            protected_indices.add(i)

    # Rest
    eff_homo = max(max_homo, 4)
    if crispr_mode and rest_tail:
        tail_len = len(rest_tail)
        body_len = max(0, rest_len - tail_len)
        while True:
            if body_len > 1:
                body_pref = "".join(random.choices(["A", "C", "G", "T"], k=body_len - 1))
            else:
                body_pref = ""
            if body_len >= 1:
                # If rest_tail starts with T, ensure base before rest_tail is NOT T
                body_last = random.choice(["A", "C", "G"]) if rest_tail.startswith("T") else random.choice(["A", "C", "G", "T"])
                cand_rest = body_pref + body_last + rest_tail
            else:
                cand_rest = rest_tail
            if get_max_homopolymer(cand_rest)[0] < eff_homo:
                rest = cand_rest
                break
    else:
        rest = "".join(random.choices(["A", "C", "G", "T"], k=rest_len))
        while get_max_homopolymer(rest)[0] >= max_homo:
            rest = "".join(random.choices(["A", "C", "G", "T"], k=rest_len))

    backbone = spacer1 + rca_binding + spacer2 + probe_binding + spacer3 + standalone_re_seq + rest

    # Protect rest_tail indices in CRISPR mode so optimizer does not mutate PAM
    if crispr_mode and rest_tail:
        tail_start = len(backbone) - len(rest_tail)
        for i in range(tail_start, len(backbone)):
            protected_indices.add(i)

    return backbone, protected_indices


def parse_segmented_backbone(
    backbone_seq: str,
    spacer1_len: int = 4,
    rca_primer_binding_len: int = 20,
    spacer2_len: int = 4,
    probe_binding_len: int = 20,
    spacer3_len: int = 0,
    rest_len: int = 4,
    restriction_site_len: int = 0,
    min_tm: float = 40.0,
    re_mode: str = "standalone",
) -> dict:
    """Parses an assembled backbone sequence into its functional sub-segments
    and derives the corresponding 5'-Biotin RCA Primer and 5'-FAM Detection Probe.
    """
    bb = str(backbone_seq).strip().upper()
    
    actual_re_len = restriction_site_len if re_mode == "standalone" else 0

    s1_end = spacer1_len
    rca_end = s1_end + rca_primer_binding_len
    s2_end = rca_end + spacer2_len
    pr_end = s2_end + probe_binding_len
    s3_end = pr_end + spacer3_len
    re_end = s3_end + actual_re_len

    spacer1 = bb[:s1_end]
    rca_binding = bb[s1_end:rca_end]
    spacer2 = bb[rca_end:s2_end]
    probe_binding = bb[s2_end:pr_end]
    spacer3 = bb[pr_end:s3_end] if spacer3_len > 0 else ""
    re_seq = bb[s3_end:re_end] if actual_re_len > 0 else ""
    rest = bb[re_end:]

    # RCA Primer binds directly to RCA_Primer_Binding_Region on circular PLP
    rca_primer_seq = str(Seq(rca_binding).reverse_complement())
    rca_gc, rca_tm = compute_metrics(Seq(rca_primer_seq))
    rca_homo_len, rca_homo_seq = get_max_homopolymer(rca_primer_seq)
    rca_hp_found, rca_hp_tm, _ = check_hairpin(rca_primer_seq)

    # Detection Probe binds to newly synthesized concatemer (Reverse_Complement(PLP))
    probe_seq = probe_binding
    probe_gc, probe_tm = compute_metrics(Seq(probe_seq))
    probe_homo_len, probe_homo_seq = get_max_homopolymer(probe_seq)
    probe_hp_found, probe_hp_tm, _ = check_hairpin(probe_seq)

    rca_tm_passes = rca_tm >= min_tm
    probe_tm_passes = probe_tm >= min_tm

    return {
        "Spacer1_Seq": spacer1,
        "RCA_Primer_Binding_Region": rca_binding,
        "Spacer2_Seq": spacer2,
        "Probe_Binding_Region": probe_binding,
        "Spacer3_Seq": spacer3,
        "Restriction_Site_Seq": re_seq,
        "Rest_Backbone_Seq": rest,
        "RCA_Primer_Seq_5to3": rca_primer_seq,
        "RCA_Primer_5_Mod": "Biotin",
        "RCA_Primer_Tm": rca_tm,
        "RCA_Primer_GC%": rca_gc,
        "RCA_Primer_Homopolymer": f"{rca_homo_len} ({rca_homo_seq})",
        "RCA_Primer_Hairpin_Tm": rca_hp_tm,
        "RCA_Primer_Valid_Tm": rca_tm_passes,
        "Detection_Probe_Seq_5to3": probe_seq,
        "Detection_Probe_5_Mod": "6-FAM",
        "Detection_Probe_Tm": probe_tm,
        "Detection_Probe_GC%": probe_gc,
        "Detection_Probe_Homopolymer": f"{probe_homo_len} ({probe_homo_seq})",
        "Detection_Probe_Hairpin_Tm": probe_hp_tm,
        "Detection_Probe_Valid_Tm": probe_tm_passes,
    }
