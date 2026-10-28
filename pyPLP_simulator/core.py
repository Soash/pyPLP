#!/usr/bin/env python3
"""
pyPLP_simulator.core: Core simulation and biophysical validation engine for
in silico padlock probe (PLP) hybridization, circularization, RCA priming,
concatemer generation, and detection probe binding.
"""

import math
import re
from typing import Dict, List, Optional, Tuple, Any
from Bio.Seq import Seq

from pyPLP.thermodynamics.tm import calculate_tm
from pyPLP.thermodynamics.metrics import compute_metrics
from pyPLP.thermodynamics.hairpin import check_hairpin
from pyPLP.thermodynamics.homopolymer import get_max_homopolymer
from pyPLP.arm_design.dimer_eval import eval_arm_dimer_vienna
from pyPLP.plp_generation.segmented_backbone import resolve_restriction_site, COMMON_RESTRICTION_ENZYMES


def compute_dimer_dg(seq1: str, seq2: str, temperature: float = 37.0) -> float:
    """Calculates heterodimer / homodimer Delta G (kcal/mol) using ViennaRNA."""
    try:
        _, dg = eval_arm_dimer_vienna(seq1, seq2, temperature=temperature)
        return round(float(dg), 2)
    except Exception:
        return 0.0


def clean_sequence(seq: str, name: str = "Sequence") -> str:
    """Sanitizes nucleic acid sequence string, strips whitespace, and converts to uppercase."""
    if seq is None:
        raise ValueError(f"Error: {name} is None.")
    s = re.sub(r"\s+", "", str(seq)).strip().upper()
    clean = re.sub(r"[^ACGTU]", "", s)
    if not clean:
        raise ValueError(f"Error: {name} contains no valid nucleotides.")
    return clean


def to_dna(seq: str) -> str:
    """Converts RNA sequence (U) to DNA representation (T) for unified processing."""
    return seq.replace("U", "T")


def reverse_complement(seq: str) -> str:
    """Computes reverse complement of DNA/RNA sequence."""
    dna_seq = to_dna(seq)
    return str(Seq(dna_seq).reverse_complement())


def calculate_duplex_dg_37(seq: str) -> float:
    """
    Approximates duplex binding Delta G (kcal/mol) at 37°C using nearest-neighbor stacking.
    """
    clean = to_dna(seq)
    if len(clean) < 2:
        return 0.0
    nn_dg37 = {
        "AA": -1.00, "TT": -1.00, "AT": -0.88, "TA": -0.58,
        "CA": -1.45, "TG": -1.45, "GT": -1.44, "AC": -1.44,
        "CT": -1.28, "AG": -1.28, "GA": -1.30, "TC": -1.30,
        "CG": -2.17, "GC": -2.24, "GG": -1.84, "CC": -1.84,
    }
    total_dg = 0.98 + 0.98  # Initiation penalty
    if clean[0] in "AT":
        total_dg += 0.05
    if clean[-1] in "AT":
        total_dg += 0.05
    for i in range(len(clean) - 1):
        dinuc = clean[i:i+2]
        total_dg += nn_dg37.get(dinuc, -1.2)
    return round(total_dg, 2)


class TargetBindingResult:
    def __init__(self):
        self.left_arm_seq: str = ""
        self.right_arm_seq: str = ""
        self.backbone_seq: str = ""
        self.left_arm_len: int = 0
        self.right_arm_len: int = 0
        self.backbone_len: int = 0
        
        self.target_left_start: int = -1
        self.target_left_end: int = -1
        self.target_right_start: int = -1
        self.target_right_end: int = -1
        self.gap_bp: int = 0
        self.is_flush_nick: bool = False
        self.junction_nick_pos: int = -1
        
        self.left_arm_tm: float = 0.0
        self.right_arm_tm: float = 0.0
        self.tm_diff: float = 0.0
        self.left_arm_dg: float = 0.0
        self.right_arm_dg: float = 0.0
        self.aggregate_tm: float = 0.0
        
        self.left_mismatches: List[Tuple[int, str, str]] = []  # (pos_in_arm, plp_base, target_comp_base)
        self.right_mismatches: List[Tuple[int, str, str]] = []
        self.nick_minus1_match: bool = False
        self.nick_plus1_match: bool = False
        
        self.alignment_display: List[str] = []
        self.verdict: str = "PASS"
        self.notes: List[str] = []


def find_best_arm_mapping(
    target_dna: str,
    plp_dna: str,
    min_arm_len: int = 10,
    max_arm_len: int = 30,
) -> TargetBindingResult:
    """
    Deconstructs PLP into Left Arm, Central Backbone, and Right Arm by identifying
    their exact binding regions on the target sequence.
    """
    res = TargetBindingResult()
    len_plp = len(plp_dna)
    len_tgt = len(target_dna)
    
    if len_plp < (min_arm_len * 2):
        raise ValueError(f"PLP length ({len_plp} nt) is too short to contain two arms of >= {min_arm_len} nt.")
    
    best_score = -999999
    best_candidate = None
    limit_max_arm = min(max_arm_len, len_plp // 2 - 5)
    
    # Target is 5'->3'. Probe arms hybridize antiparallel.
    # Left Arm (5'->3') hybridizes to Target_Left (5'->3'): revcomp(Left Arm) == Target_Left
    # Right Arm (5'->3') hybridizes to Target_Right (5'->3'): revcomp(Right Arm) == Target_Right
    for l_len in range(min_arm_len, limit_max_arm + 1):
        left_cand = plp_dna[:l_len]
        left_target_comp = reverse_complement(left_cand)
        
        for match_l in re.finditer(re.escape(left_target_comp), target_dna):
            pos_l_start = match_l.start()
            pos_l_end = match_l.end()
            
            for r_len in range(min_arm_len, limit_max_arm + 1):
                right_cand = plp_dna[len_plp - r_len:]
                right_target_comp = reverse_complement(right_cand)
                
                pos_r_start = pos_l_end
                pos_r_end = pos_r_start + r_len
                
                if pos_r_end <= len(target_dna) and target_dna[pos_r_start:pos_r_end] == right_target_comp:
                    score = (l_len + r_len) * 10
                    if score > best_score:
                        best_score = score
                        best_candidate = {
                            "l_len": l_len,
                            "r_len": r_len,
                            "pos_l_start": pos_l_start,
                            "pos_l_end": pos_l_end,
                            "pos_r_start": pos_r_start,
                            "pos_r_end": pos_r_end,
                            "gap": 0,
                            "exact": True,
                        }
    
    if not best_candidate:
        best_l = None
        for l_len in range(15, 11, -1):
            left_cand = plp_dna[:l_len]
            t_comp = reverse_complement(left_cand)
            pos = target_dna.find(t_comp)
            if pos != -1:
                best_l = (l_len, pos, pos + l_len)
                break
        
        best_r = None
        for r_len in range(15, 11, -1):
            right_cand = plp_dna[len_plp - r_len:]
            t_comp = reverse_complement(right_cand)
            pos = target_dna.find(t_comp)
            if pos != -1:
                best_r = (r_len, pos, pos + r_len)
                break
        
        if best_l and best_r:
            gap = best_r[1] - best_l[2]
            best_candidate = {
                "l_len": best_l[0],
                "r_len": best_r[0],
                "pos_l_start": best_l[1],
                "pos_l_end": best_l[2],
                "pos_r_start": best_r[1],
                "pos_r_end": best_r[2],
                "gap": gap,
                "exact": False,
            }
        else:
            def_arm = 15 if len_plp >= 40 else len_plp // 3
            best_candidate = {
                "l_len": def_arm,
                "r_len": def_arm,
                "pos_l_start": 0,
                "pos_l_end": def_arm,
                "pos_r_start": def_arm,
                "pos_r_end": def_arm * 2,
                "gap": 0,
                "exact": False,
            }

    l_len = best_candidate["l_len"]
    r_len = best_candidate["r_len"]
    res.left_arm_seq = plp_dna[:l_len]
    res.right_arm_seq = plp_dna[len_plp - r_len:]
    res.backbone_seq = plp_dna[l_len : len_plp - r_len]
    
    res.left_arm_len = l_len
    res.right_arm_len = r_len
    res.backbone_len = len(res.backbone_seq)
    
    res.target_left_start = best_candidate["pos_l_start"]
    res.target_left_end = best_candidate["pos_l_end"]
    res.target_right_start = best_candidate["pos_r_start"]
    res.target_right_end = best_candidate["pos_r_end"]
    res.gap_bp = best_candidate["gap"]
    res.is_flush_nick = (res.gap_bp == 0)
    res.junction_nick_pos = res.target_left_end
    
    _, res.left_arm_tm = compute_metrics(Seq(res.left_arm_seq))
    _, res.right_arm_tm = compute_metrics(Seq(res.right_arm_seq))
    res.tm_diff = round(abs(res.right_arm_tm - res.left_arm_tm), 2)
    res.left_arm_dg = calculate_duplex_dg_37(res.left_arm_seq)
    res.right_arm_dg = calculate_duplex_dg_37(res.right_arm_seq)
    
    full_target_hybrid = res.left_arm_seq + res.right_arm_seq
    _, res.aggregate_tm = compute_metrics(Seq(full_target_hybrid))
    
    # Evaluate Mismatches
    if 0 <= res.target_right_start < len(target_dna) and res.target_right_end <= len(target_dna):
        t_sub_r = target_dna[res.target_right_start:res.target_right_end]
        exp_r = reverse_complement(t_sub_r)
        for i, (p_char, e_char) in enumerate(zip(res.right_arm_seq, exp_r)):
            if p_char != e_char:
                res.right_mismatches.append((i + 1, p_char, e_char))
        if res.right_arm_seq and exp_r:
            res.nick_minus1_match = (res.right_arm_seq[-1] == exp_r[-1])
    else:
        res.nick_minus1_match = True
        
    if 0 <= res.target_left_start < len(target_dna) and res.target_left_end <= len(target_dna):
        t_sub_l = target_dna[res.target_left_start:res.target_left_end]
        exp_l = reverse_complement(t_sub_l)
        for i, (p_char, e_char) in enumerate(zip(res.left_arm_seq, exp_l)):
            if p_char != e_char:
                res.left_mismatches.append((i + 1, p_char, e_char))
        if res.left_arm_seq and exp_l:
            res.nick_plus1_match = (res.left_arm_seq[0] == exp_l[0])
    else:
        res.nick_plus1_match = True

    if not res.is_flush_nick:
        res.verdict = "FAIL"
        res.notes.append(f"Gap of {res.gap_bp} bp detected between arm binding sites! Padlock circularization cannot proceed without a flush nick.")
    elif res.left_mismatches or res.right_mismatches:
        if not res.nick_minus1_match:
            res.verdict = "FAIL"
            res.notes.append("Critical mismatch at position -1 (3' ligation junction terminus)! DNA ligase cannot circularize this probe.")
        elif not res.nick_plus1_match:
            res.verdict = "WARNING"
            res.notes.append("Mismatch at position +1 (5' junction terminus) - ligation efficiency may be reduced.")
        else:
            res.verdict = "WARNING"
            res.notes.append(f"Internal mismatch(es) detected: Left={len(res.left_mismatches)}, Right={len(res.right_mismatches)}.")
    else:
        res.verdict = "PASS"
        res.notes.append("Perfect 0-bp flush nick alignment with 100% sequence complementarity across both arms.")

    return res


class BackboneAnnotationResult:
    def __init__(self):
        self.spacer1_seq: str = ""
        self.rca_primer_binding_seq: str = ""
        self.spacer2_seq: str = ""
        self.probe_binding_seq: str = ""
        self.restriction_site_seq: str = ""
        self.restriction_enzyme_name: str = ""
        self.rest_seq: str = ""
        
        self.rca_primer_match_found: bool = False
        self.rca_primer_start: int = -1
        self.rca_primer_end: int = -1
        
        self.probe_match_found: bool = False
        self.probe_start: int = -1
        self.probe_end: int = -1
        self.probe_binding_mode: str = "CONCATEMER"
        
        self.notes: List[str] = []


def annotate_backbone_regions(
    backbone_seq: str,
    rca_primer_seq: Optional[str] = None,
    probe_seq: Optional[str] = None,
) -> BackboneAnnotationResult:
    """
    Scans the PLP central backbone to locate the RCA Primer Binding Site,
    Detection Probe Binding Site, Restriction Cleavage Site, and Spacer Regions.
    """
    res = BackboneAnnotationResult()
    bb = to_dna(clean_sequence(backbone_seq, "Backbone"))
    bb_len = len(bb)
    
    if rca_primer_seq:
        clean_primer = to_dna(clean_sequence(rca_primer_seq, "RCA Primer"))
        target_rca_site = reverse_complement(clean_primer)
        
        idx = bb.find(target_rca_site)
        if idx != -1:
            res.rca_primer_match_found = True
            res.rca_primer_start = idx
            res.rca_primer_end = idx + len(target_rca_site)
            res.rca_primer_binding_seq = target_rca_site
            res.notes.append(f"RCA Primer binding site identified at backbone nt {idx+1}-{idx+len(target_rca_site)}.")
        else:
            for p_len in range(len(clean_primer), 10, -1):
                sub_primer = clean_primer[:p_len]
                sub_site = reverse_complement(sub_primer)
                idx = bb.find(sub_site)
                if idx != -1:
                    res.rca_primer_match_found = True
                    res.rca_primer_start = idx
                    res.rca_primer_end = idx + len(sub_site)
                    res.rca_primer_binding_seq = sub_site
                    res.notes.append(f"Partial RCA Primer binding site ({p_len}/{len(clean_primer)} nt) identified at nt {idx+1}-{idx+len(sub_site)}.")
                    break
    
    if probe_seq:
        clean_pr = to_dna(clean_sequence(probe_seq, "Detection Probe"))
        idx = bb.find(clean_pr)
        if idx != -1:
            res.probe_match_found = True
            res.probe_start = idx
            res.probe_end = idx + len(clean_pr)
            res.probe_binding_seq = clean_pr
            res.probe_binding_mode = "CONCATEMER"
            res.notes.append(f"Concatemer Detection Probe cassette found at backbone nt {idx+1}-{idx+len(clean_pr)} (Identical motif; hybridizes to concatemer).")
        else:
            pr_rc = reverse_complement(clean_pr)
            idx_rc = bb.find(pr_rc)
            if idx_rc != -1:
                res.probe_match_found = True
                res.probe_start = idx_rc
                res.probe_end = idx_rc + len(pr_rc)
                res.probe_binding_seq = pr_rc
                res.probe_binding_mode = "DIRECT_PLP"
                res.notes.append(f"Direct PLP Detection Probe site found at backbone nt {idx_rc+1}-{idx_rc+len(pr_rc)} (Direct hybrid to PLP circle).")
    
    for enz_name, enz_seq in sorted(COMMON_RESTRICTION_ENZYMES.items(), key=lambda x: -len(x[1])):
        idx_re = bb.find(enz_seq)
        if idx_re != -1:
            res.restriction_enzyme_name = enz_name
            res.restriction_site_seq = enz_seq
            res.notes.append(f"Restriction endonuclease cleavage motif '{enz_name}' ({enz_seq}) found at nt {idx_re+1}-{idx_re+len(enz_seq)}.")
            break
            
    if res.rca_primer_match_found and res.probe_match_found:
        if res.rca_primer_start < res.probe_start:
            res.spacer1_seq = bb[:res.rca_primer_start]
            res.spacer2_seq = bb[res.rca_primer_end : res.probe_start]
            res.rest_seq = bb[res.probe_end:]
        else:
            res.spacer1_seq = bb[:res.probe_start]
            res.spacer2_seq = bb[res.probe_end : res.rca_primer_start]
            res.rest_seq = bb[res.rca_primer_end:]
    elif res.rca_primer_match_found:
        res.spacer1_seq = bb[:res.rca_primer_start]
        res.rest_seq = bb[res.rca_primer_end:]
    else:
        s1 = 4
        rca_len = 15 if bb_len >= 40 else bb_len // 3
        s2 = 4
        pr_len = 15 if bb_len >= 40 else bb_len // 3
        res.spacer1_seq = bb[:s1]
        res.rca_primer_binding_seq = bb[s1 : s1 + rca_len]
        res.spacer2_seq = bb[s1 + rca_len : s1 + rca_len + s2]
        res.probe_binding_seq = bb[s1 + rca_len + s2 : s1 + rca_len + s2 + pr_len]
        res.rest_seq = bb[s1 + rca_len + s2 + pr_len :]
        
    return res


class RCAPrimingResult:
    def __init__(self):
        self.primer_seq: str = ""
        self.primer_len: int = 0
        self.primer_tm: float = 0.0
        self.primer_gc: float = 0.0
        self.primer_dg: float = 0.0
        self.binding_site_on_plp: str = ""
        self.priming_direction: str = "5' -> 3' Rolling Circle Synthesis"
        self.extension_strand: str = "Concatemer (Reverse Complement of PLP)"
        self.verdict: str = "PASS"
        self.notes: List[str] = []


def simulate_rca_priming(
    plp_seq: str,
    rca_primer_seq: str,
    reaction_temp: float = 37.0,
) -> RCAPrimingResult:
    """Simulates annealing of the RCA Primer onto the circularized padlock probe."""
    res = RCAPrimingResult()
    p_clean = to_dna(clean_sequence(rca_primer_seq, "RCA Primer"))
    plp_clean = to_dna(clean_sequence(plp_seq, "Padlock Probe"))
    
    res.primer_seq = p_clean
    res.primer_len = len(p_clean)
    res.primer_gc, res.primer_tm = compute_metrics(Seq(p_clean))
    res.primer_dg = calculate_duplex_dg_37(p_clean)
    
    expected_site = reverse_complement(p_clean)
    idx = plp_clean.find(expected_site)
    
    if idx != -1:
        res.binding_site_on_plp = expected_site
        res.notes.append(f"RCA Primer hybridizes perfectly to circular PLP coordinates {idx+1}-{idx+len(expected_site)}.")
        if res.primer_tm < reaction_temp:
            res.verdict = "WARNING"
            res.notes.append(f"RCA Primer Tm ({res.primer_tm}°C) is below reaction temperature ({reaction_temp}°C).")
        else:
            res.verdict = "PASS"
            res.notes.append(f"Primer Tm ({res.primer_tm}°C) is comfortably above reaction temperature ({reaction_temp}°C). Efficient strand displacement initiation expected.")
    else:
        res.verdict = "WARNING"
        res.notes.append("Exact full-length primer binding cassette not detected in PLP sequence. Check primer reverse-complement orientation.")
        
    return res


class ConcatemerResult:
    def __init__(self):
        self.unit_repeat_seq: str = ""
        self.unit_repeat_len: int = 0
        self.concatemer_3x_seq: str = ""
        self.repeating_probe_sites_count: int = 0
        self.probe_binding_positions: List[int] = []
        self.probe_tm: float = 0.0
        self.probe_gc: float = 0.0
        self.probe_dg: float = 0.0
        self.verdict: str = "PASS"
        self.notes: List[str] = []


def generate_rca_concatemer(
    plp_seq: str,
    probe_seq: Optional[str] = None,
    num_repeats: int = 3,
    reaction_temp: float = 37.0,
) -> ConcatemerResult:
    """
    Generates rolling circle amplification (RCA) concatemer sequence:
    Each turn of Phi29 polymerase generates 1 unit of Reverse_Complement(PLP).
    """
    res = ConcatemerResult()
    plp_clean = to_dna(clean_sequence(plp_seq, "Padlock Probe"))
    
    unit = reverse_complement(plp_clean)
    res.unit_repeat_seq = unit
    res.unit_repeat_len = len(unit)
    res.concatemer_3x_seq = unit * num_repeats
    
    if probe_seq:
        pr_clean = to_dna(clean_sequence(probe_seq, "Detection Probe"))
        res.probe_gc, res.probe_tm = compute_metrics(Seq(pr_clean))
        res.probe_dg = calculate_duplex_dg_37(pr_clean)
        
        target_motif_on_concatemer = reverse_complement(pr_clean)
        matches = [m.start() for m in re.finditer(re.escape(target_motif_on_concatemer), res.concatemer_3x_seq)]
        res.probe_binding_positions = matches
        res.repeating_probe_sites_count = len(matches)
        
        if res.repeating_probe_sites_count >= num_repeats:
            res.verdict = "PASS"
            res.notes.append(f"Detection probe binds exactly {res.repeating_probe_sites_count // num_repeats} time(s) per single RCA monomer unit.")
            res.notes.append(f"In 3 repeats, {res.repeating_probe_sites_count} tandem reporter fluorophore binding sites are generated.")
        elif res.repeating_probe_sites_count > 0:
            res.verdict = "PASS"
            res.notes.append(f"{res.repeating_probe_sites_count} probe hybridization sites detected across 3 repeat monomers.")
        else:
            direct_matches = [m.start() for m in re.finditer(re.escape(pr_clean), res.concatemer_3x_seq)]
            if direct_matches:
                res.verdict = "PASS"
                res.probe_binding_positions = direct_matches
                res.repeating_probe_sites_count = len(direct_matches)
                res.notes.append(f"Probe sequence matches directly on concatemer ({res.repeating_probe_sites_count} sites).")
            else:
                res.verdict = "WARNING"
                res.notes.append("Detection probe target sequence not detected within concatemer repeat units.")
                
        if res.probe_tm < reaction_temp:
            res.notes.append(f"Notice: Probe Tm ({res.probe_tm}°C) is below reaction temperature ({reaction_temp}°C).")
            
    return res


class QCMatrixResult:
    def __init__(self):
        self.plp_hairpin_tm: float = 0.0
        self.plp_hairpin_found: str = "No"
        self.primer_hairpin_tm: float = 0.0
        self.primer_hairpin_found: str = "No"
        self.probe_hairpin_tm: float = 0.0
        self.probe_hairpin_found: str = "No"
        
        self.plp_self_dimer_dg: float = 0.0
        self.primer_self_dimer_dg: float = 0.0
        self.probe_self_dimer_dg: float = 0.0
        
        self.primer_probe_dimer_dg: float = 0.0
        self.plp_primer_dimer_dg: float = 0.0
        self.plp_probe_dimer_dg: float = 0.0
        
        self.verdict: str = "PASS"
        self.warnings: List[str] = []


def run_cross_reactivity_qc(
    plp_seq: str,
    rca_primer_seq: Optional[str] = None,
    probe_seq: Optional[str] = None,
    threshold_dg: float = -3.0,
) -> QCMatrixResult:
    """Evaluates unimolecular hairpins, homodimers, and pairwise heterodimerization between all oligos."""
    res = QCMatrixResult()
    plp_clean = to_dna(clean_sequence(plp_seq, "PLP"))
    
    res.plp_hairpin_found, res.plp_hairpin_tm, _ = check_hairpin(plp_clean)
    res.plp_self_dimer_dg = compute_dimer_dg(plp_clean, plp_clean)
    
    p_clean = to_dna(clean_sequence(rca_primer_seq, "RCA Primer")) if rca_primer_seq else ""
    pr_clean = to_dna(clean_sequence(probe_seq, "Probe")) if probe_seq else ""
    
    if p_clean:
        res.primer_hairpin_found, res.primer_hairpin_tm, _ = check_hairpin(p_clean)
        res.primer_self_dimer_dg = compute_dimer_dg(p_clean, p_clean)
        res.plp_primer_dimer_dg = compute_dimer_dg(plp_clean, p_clean)
        
    if pr_clean:
        res.probe_hairpin_found, res.probe_hairpin_tm, _ = check_hairpin(pr_clean)
        res.probe_self_dimer_dg = compute_dimer_dg(pr_clean, pr_clean)
        res.plp_probe_dimer_dg = compute_dimer_dg(plp_clean, pr_clean)
        
    if p_clean and pr_clean:
        res.primer_probe_dimer_dg = compute_dimer_dg(p_clean, pr_clean)
        if res.primer_probe_dimer_dg <= threshold_dg:
            res.warnings.append(f"Significant cross-dimer between RCA Primer and Detection Probe (dG = {res.primer_probe_dimer_dg} kcal/mol). Risk of primer-probe titration.")
            res.verdict = "WARNING"

    if res.primer_self_dimer_dg <= threshold_dg:
        res.warnings.append(f"RCA Primer self-dimer detected (dG = {res.primer_self_dimer_dg} kcal/mol).")
        res.verdict = "WARNING"
        
    if res.probe_self_dimer_dg <= threshold_dg:
        res.warnings.append(f"Detection Probe self-dimer detected (dG = {res.probe_self_dimer_dg} kcal/mol).")
        res.verdict = "WARNING"
        
    return res


class SimulationReport:
    """Master simulation report encapsulating all 6 reaction stages."""
    def __init__(
        self,
        target_seq: str,
        plp_seq: str,
        rca_primer_seq: Optional[str] = None,
        probe_seq: Optional[str] = None,
        reaction_temp: float = 37.0,
        is_rna: bool = False,
    ):
        self.raw_target: str = target_seq
        self.target_dna: str = to_dna(clean_sequence(target_seq, "Target"))
        self.plp_dna: str = to_dna(clean_sequence(plp_seq, "Padlock Probe"))
        self.rca_primer_dna: str = to_dna(clean_sequence(rca_primer_seq, "RCA Primer")) if rca_primer_seq else ""
        self.probe_dna: str = to_dna(clean_sequence(probe_seq, "Detection Probe")) if probe_seq else ""
        self.reaction_temp: float = reaction_temp
        self.is_rna: bool = is_rna or ("U" in target_seq.upper())
        
        self.target_binding: TargetBindingResult = find_best_arm_mapping(self.target_dna, self.plp_dna)
        self.backbone_annotation: BackboneAnnotationResult = annotate_backbone_regions(
            self.target_binding.backbone_seq,
            self.rca_primer_dna,
            self.probe_dna,
        )
        self.rca_priming: Optional[RCAPrimingResult] = None
        if self.rca_primer_dna:
            self.rca_priming = simulate_rca_priming(self.plp_dna, self.rca_primer_dna, reaction_temp)
            
        self.concatemer: ConcatemerResult = generate_rca_concatemer(
            self.plp_dna,
            self.probe_dna if self.probe_dna else None,
            num_repeats=3,
            reaction_temp=reaction_temp,
        )
        
        self.qc_matrix: QCMatrixResult = run_cross_reactivity_qc(
            self.plp_dna,
            self.rca_primer_dna if self.rca_primer_dna else None,
            self.probe_dna if self.probe_dna else None,
        )
        
        self.overall_verdict: str = "PASS"
        if self.target_binding.verdict == "FAIL":
            self.overall_verdict = "FAIL"
        elif any(v.verdict == "WARNING" for v in [self.target_binding, self.qc_matrix] if v):
            self.overall_verdict = "WARNING"
