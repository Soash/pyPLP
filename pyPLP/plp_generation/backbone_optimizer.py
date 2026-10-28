from Bio.Seq import Seq

from pyPLP.plp_generation.random_backbone import generate_random_backbone
from pyPLP.plp_generation.segmented_backbone import (
    build_initial_segmented_backbone,
    parse_segmented_backbone,
    resolve_restriction_site,
)
from pyPLP.plp_generation.stem_dissolution import break_stem_mutations_aggressive
from pyPLP.plp_generation.structure_eval import eval_plp_structure
from pyPLP.thermodynamics.hairpin import check_hairpin
from pyPLP.thermodynamics.homopolymer import get_max_homopolymer
from pyPLP.thermodynamics.metrics import compute_metrics


def optimize_backbone_candidate(
    initial_bb: str,
    left_arm: str,
    right_arm: str,
    max_homo: int = 4,
    max_hairpin_tm: float = 35.0,
    min_plp_mfe_dg: float = -6.0,
    max_paired_arm_bases: int = 4,
    max_paired_total_bases: int = 16,
    folding_engine: str = "RNAfold",
    max_cycles: int = 30,
    arch_cfg: dict | None = None,
    protected_bb_indices: set[int] | None = None,
) -> dict | None:
    """Stochastically optimizes a backbone sequence for a given left and right arm pair
    to guarantee high openness, minimal secondary structure, and compliant LFA oligos.
    """
    len_l = len(left_arm)
    len_r = len(right_arm)
    curr_bb = initial_bb
    cur_prot_indices = set(protected_bb_indices) if protected_bb_indices else set()

    def _regenerate_bb():
        if arch_cfg:
            return build_initial_segmented_backbone(
                spacer1_len=arch_cfg.get("spacer1_len", 4),
                rca_primer_binding_len=arch_cfg.get("rca_primer_binding_len", 20),
                spacer2_len=arch_cfg.get("spacer2_len", 4),
                probe_binding_len=arch_cfg.get("probe_binding_len", 20),
                spacer3_len=arch_cfg.get("spacer3_len", 0),
                rest_len=arch_cfg.get("rest_len", 4),
                restriction_site=arch_cfg.get("restriction_site", ""),
                universal_rca_primer_seq=arch_cfg.get("universal_rca_primer_seq", ""),
                universal_probe_seq=arch_cfg.get("universal_probe_seq", ""),
                min_tm=arch_cfg.get("min_rca_primer_tm", 40.0),
                max_homo=max_homo,
                max_hairpin_tm=max_hairpin_tm,
                re_mode=arch_cfg.get("re_mode", "standalone"),
                re_probe_offset=arch_cfg.get("re_probe_offset", 2),
                crispr_mode=bool(arch_cfg.get("crispr_mode", False)),
                rest_tail=arch_cfg.get("rest_tail", ""),
            )
        else:
            return generate_random_backbone(len(curr_bb)), set()

    crispr_active = bool(arch_cfg.get("crispr_mode", False)) if arch_cfg else False
    eff_homo = max(max_homo, 4) if crispr_active else max_homo

    for _ in range(max_cycles):
        full_plp = left_arm + curr_bb + right_arm

        bb_homo_len, _ = get_max_homopolymer(curr_bb)
        if bb_homo_len >= max_homo:
            curr_bb, cur_prot_indices = _regenerate_bb()
            continue

        jl_len, _ = get_max_homopolymer(left_arm[-3:] + curr_bb[:3])
        jr_len, _ = get_max_homopolymer(curr_bb[-3:] + right_arm[:3])
        if jl_len >= max_homo or jr_len >= eff_homo:
            curr_bb, cur_prot_indices = _regenerate_bb()
            continue

        bb_hp_found, bb_hp_tm, _ = check_hairpin(curr_bb)
        if bb_hp_found == "Yes" and bb_hp_tm > max_hairpin_tm:
            #print(f"Rejected hairpin: {bb_hp_tm}")
            curr_bb, cur_prot_indices = _regenerate_bb()
            continue
            curr_bb, cur_prot_indices = _regenerate_bb()
            continue

        dot_bracket, mfe_dg = eval_plp_structure(full_plp, folding_engine=folding_engine)

        paired_left = sum(1 for ch in dot_bracket[:len_l] if ch in "()")
        paired_right = sum(1 for ch in dot_bracket[-len_r:] if ch in "()")
        paired_total = sum(1 for ch in dot_bracket if ch in "()")

        is_open = (
            mfe_dg >= min_plp_mfe_dg
            and paired_left <= max_paired_arm_bases
            and paired_right <= max_paired_arm_bases
            and paired_total <= max_paired_total_bases
        )

        if is_open and arch_cfg and arch_cfg.get("screen_against"):
            from pyPLP.multiplex.cross_dimer import calc_heterodimer_dg
            for s_seq in arch_cfg["screen_against"]:
                if calc_heterodimer_dg(curr_bb, s_seq) < arch_cfg.get("max_cross_dimer_dg", -5.0):
                    is_open = False
                    curr_bb, cur_prot_indices = _regenerate_bb()
                    break
            
            if not is_open:
                continue

        if not is_open: pass
        if is_open:
            seg_data = {}
            if arch_cfg is not None:
                re_site_str = arch_cfg.get("restriction_site", "")
                re_len = len(resolve_restriction_site(re_site_str)) if re_site_str else 0
                seg_data = parse_segmented_backbone(
                    curr_bb,
                    spacer1_len=arch_cfg.get("spacer1_len", 4),
                    rca_primer_binding_len=arch_cfg.get("rca_primer_binding_len", 20),
                    spacer2_len=arch_cfg.get("spacer2_len", 4),
                    probe_binding_len=arch_cfg.get("probe_binding_len", 20),
                    spacer3_len=arch_cfg.get("spacer3_len", 0),
                    rest_len=arch_cfg.get("rest_len", 4),
                    restriction_site_len=re_len,
                    min_tm=arch_cfg.get("min_rca_primer_tm", 40.0),
                    re_mode=arch_cfg.get("re_mode", "standalone"),
                )
                if not seg_data.get("RCA_Primer_Valid_Tm", True):
                    curr_bb, cur_prot_indices = _regenerate_bb()
                    continue
                if not seg_data.get("Detection_Probe_Valid_Tm", True):
                    curr_bb, cur_prot_indices = _regenerate_bb()
                    continue

            bb_homo_len, bb_homo_seq = get_max_homopolymer(curr_bb)
            bb_gc, bb_tm = compute_metrics(Seq(curr_bb))
            full_gc, full_tm = compute_metrics(Seq(full_plp))

            res = {
                "bb": curr_bb,
                "full_plp": full_plp,
                "full_gc": full_gc,
                "full_tm": full_tm,
                "full_tm_std": full_tm,
                "full_tm_pcr": full_tm,
                "bb_gc": bb_gc,
                "bb_tm": bb_tm,
                "bb_tm_std": bb_tm,
                "bb_tm_pcr": bb_tm,
                "bb_homo": f"{bb_homo_len} ({bb_homo_seq})",
                "bb_hp_found": bb_hp_found,
                "bb_hp_tm": bb_hp_tm,
                "dot_bracket": dot_bracket,
                "mfe_dg": mfe_dg,
                "paired_left": paired_left,
                "paired_right": paired_right,
                "paired_total": paired_total,
            }
            res.update(seg_data)
            return res

        prot_plp_indices = {idx + len_l for idx in cur_prot_indices} if cur_prot_indices else None
        mutated_plp, was_mutated = break_stem_mutations_aggressive(
            full_plp, dot_bracket, len_l, len_r, max_homo=max_homo, protected_indices=prot_plp_indices
        )
        if not was_mutated:
            break

        curr_bb = mutated_plp[len_l:-len_r]

    return None
