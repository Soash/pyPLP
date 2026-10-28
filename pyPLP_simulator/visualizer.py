#!/usr/bin/env python3
"""
pyPLP_simulator.visualizer: Generates detailed terminal and ASCII diagrams for
padlock probe target hybridization, ligation nick junction, circular topology,
RCA priming, and concatemer detection probe hybridization.
"""

from Bio.Seq import Seq
from pyPLP_simulator.core import SimulationReport, reverse_complement


def generate_ascii_alignment(report: SimulationReport) -> str:
    """
    Renders nucleotide-level base-pairing ladder between Target and Padlock Probe arms.
    """
    tb = report.target_binding
    target_seq = report.target_dna
    is_rna = report.is_rna
    
    l_start = tb.target_left_start
    l_end = tb.target_left_end
    r_start = tb.target_right_start
    r_end = tb.target_right_end
    
    # Extract surrounding target context (up to 10 nt on flanks)
    flank_left = max(0, l_start - 8)
    flank_right = min(len(target_seq), r_end + 8)
    
    prefix_target = target_seq[flank_left:l_start]
    target_l_part = target_seq[l_start:l_end]
    target_r_part = target_seq[r_start:r_end]
    suffix_target = target_seq[r_end:flank_right]
    
    # Format target line (RNA shows U, DNA shows T)
    disp_prefix = prefix_target.replace("T", "U") if is_rna else prefix_target
    disp_tl = target_l_part.replace("T", "U") if is_rna else target_l_part
    disp_tr = target_r_part.replace("T", "U") if is_rna else target_r_part
    disp_suffix = suffix_target.replace("T", "U") if is_rna else suffix_target
    
    target_disp = f"Target (5'->3') : ... {disp_prefix} [ {disp_tl} | {disp_tr} ] {disp_suffix} ..."
    
    # Build match symbols
    # Target Left (5'->3') pairs antiparallel with Left Arm (3'<-5'):
    # Left Arm 3' is on left, 5'-PO4 is at nick on right
    left_arm_disp = tb.left_arm_seq[::-1]
    match_l = []
    for i in range(len(target_l_part)):
        p_c = left_arm_disp[i] if i < len(left_arm_disp) else "N"
        t_c = target_l_part[i].replace("U", "T")
        match_l.append("|" if Seq(p_c).complement() == Seq(t_c) else "X")
    match_l_str = "".join(match_l)
    
    # Target Right (5'->3') pairs antiparallel with Right Arm (3'<-5'):
    # Right Arm 3'-OH is at nick on left, 5' is on right
    right_arm_disp = tb.right_arm_seq[::-1]
    match_r = []
    for i in range(len(target_r_part)):
        p_c = right_arm_disp[i] if i < len(right_arm_disp) else "N"
        t_c = target_r_part[i].replace("U", "T")
        match_r.append("|" if Seq(p_c).complement() == Seq(t_c) else "X")
    match_r_str = "".join(match_r)
    
    indent = len("Target (5'->3') : ... ") + len(disp_prefix) + 2
    pad = " " * indent
    ladder = f"{pad}  {match_l_str} ▼ {match_r_str}"
    
    # Probe line
    probe_line = f"{pad}[ 3'-{left_arm_disp}-5'-PO4 ] [ HO-3'-{right_arm_disp}-5' ]"
    
    # Nick marker line
    nick_marker = f"{pad}{' '*len(left_arm_disp)}       ▲ LIGATION NICK"
    
    lines = [
        "=" * 90,
        "STAGE 1: TARGET HYBRIDIZATION & LIGATION NICK BASE-PAIRING LADDER",
        "=" * 90,
        f"Target Type         : {'RNA (Ribonucleic Acid)' if is_rna else 'DNA (Deoxyribonucleic Acid)'}",
        f"Ligation Junction   : {'0-bp PERFECT FLUSH NICK (Ligatable)' if tb.is_flush_nick else f'FAILED ({tb.gap_bp} bp gap)'}",
        f"Nick Status (-1/+1) : 3' End (-1)={'MATCH [Ready to ligate]' if tb.nick_minus1_match else 'MISMATCH [Incompetent]'}, 5' End (+1)={'MATCH' if tb.nick_plus1_match else 'MISMATCH'}",
        f"Left Arm Tm / dG    : {tb.left_arm_tm:.1f}°C  |  {tb.left_arm_dg:.2f} kcal/mol ({tb.left_arm_len} nt)",
        f"Right Arm Tm / dG   : {tb.right_arm_tm:.1f}°C  |  {tb.right_arm_dg:.2f} kcal/mol ({tb.right_arm_len} nt)",
        f"Arm Tm Difference   : ΔTm = {tb.tm_diff:.1f}°C ({'Balanced (≤ 5°C)' if tb.tm_diff <= 5.0 else 'Asymmetric (> 5°C)'})",
        f"Total Hybrid Tm     : {tb.aggregate_tm:.1f}°C (Bivalent Duplex Footprint: {tb.left_arm_len + tb.right_arm_len} bp)",
        "-" * 90,
        target_disp,
        ladder,
        probe_line,
        nick_marker,
        "-" * 90,
    ]
    return "\n".join(lines)


def generate_circular_topology_ascii(report: SimulationReport) -> str:
    """
    Renders an ASCII visualization of the circularized padlock probe topology.
    """
    tb = report.target_binding
    bb = report.backbone_annotation
    
    lines = [
        "=" * 90,
        "STAGE 2: CIRCULAR PADLOCK PROBE TOPOLOGY & FUNCTIONAL BACKBONE ANNOTATION",
        "=" * 90,
        f"Full Padlock Length : {len(report.plp_dna)} nt",
        f"Left Hybridizing Arm: {tb.left_arm_seq} ({tb.left_arm_len} nt)",
        f"Right Hybridizing Arm: {tb.right_arm_seq} ({tb.right_arm_len} nt)",
        f"Central Connector BB : {tb.backbone_seq} ({tb.backbone_len} nt)",
        "-" * 90,
        "                      [ LIGATION PHOSPHODIESTER BOND ]                      ",
        "                                     ▼                                      ",
        f"              5'-PO4 [ Left Arm ]─────────[ Right Arm ] 3'-OH              ",
        "                       ({0:>2} nt)                 ({1:>2} nt)                ".format(tb.left_arm_len, tb.right_arm_len),
        "                         │                               │                  ",
        "                   [ Spacer 1 ]                    [ Rest / Cleave ]        ",
        "                     ({0:>2} nt)                         ({1:>2} nt)          ".format(len(bb.spacer1_seq), len(bb.rest_seq)),
        "                         │                               │                  ",
        "               [ RCA Primer Site ]              [ FAM Probe Site ]          ",
        "                     ({0:>2} nt)                         ({1:>2} nt)          ".format(len(bb.rca_primer_binding_seq), len(bb.probe_binding_seq)),
        "                         │                               │                  ",
        "                         └─────────[ Spacer 2 ]──────────┘                  ",
        "                                     ({0:>2} nt)                            ".format(len(bb.spacer2_seq)),
        "-" * 90,
        f"  • Spacer 1 Region     : {bb.spacer1_seq or 'None'} ({len(bb.spacer1_seq)} nt)",
        f"  • RCA Primer Binding  : {bb.rca_primer_binding_seq or 'None'} ({len(bb.rca_primer_binding_seq)} nt)",
        f"  • Spacer 2 Region     : {bb.spacer2_seq or 'None'} ({len(bb.spacer2_seq)} nt)",
        f"  • FAM Probe Binding   : {bb.probe_binding_seq or 'None'} ({len(bb.probe_binding_seq)} nt) [{bb.probe_binding_mode}]",
        f"  • Restriction Motif   : {bb.restriction_enzyme_name or 'None'} ({bb.restriction_site_seq or 'N/A'})",
        f"  • Rest Backbone Seq   : {bb.rest_seq or 'None'} ({len(bb.rest_seq)} nt)",
        "-" * 90,
    ]
    return "\n".join(lines)


def generate_rca_and_detection_ascii(report: SimulationReport) -> str:
    """
    Renders RCA Priming and Concatemer Detection Probe hybridization diagram.
    """
    rp = report.rca_priming
    cat = report.concatemer
    lines = [
        "=" * 90,
        "STAGE 3 & 4: RCA PRIMING & CONCATEMER DETECTION PROBE HYBRIDIZATION",
        "=" * 90,
    ]
    
    if rp:
        lines.extend([
            f"RCA Primer Sequence : 5'-{rp.primer_seq}-3' ({rp.primer_len} nt)",
            f"Primer Tm / GC%     : {rp.primer_tm:.1f}°C  |  {rp.primer_gc:.1f}% (Binding dG: {rp.primer_dg:.2f} kcal/mol)",
            f"Priming Direction   : {rp.priming_direction}",
            f"Extension Template  : Circular Single-Stranded Padlock Probe (Rolling Circle)",
            "-" * 90,
            "Primer Docking Alignment:",
            f"  Circular PLP Target : 3' ... {rp.binding_site_on_plp} ... 5'",
            f"                               {'|'*rp.primer_len}",
            f"  RCA Primer          : 5'-[Biotin]-{rp.primer_seq}-OH 3' ───► (Polymerase Extension Direction)",
            "-" * 90,
        ])
    else:
        lines.extend(["RCA Primer          : Not specified (Skipping Primer Alignment)", "-" * 90])
        
    lines.extend([
        f"Single Monomer Unit : {cat.unit_repeat_len} nt (Reverse Complement of Circular PLP)",
        f"Tandem Repeat Model : 3 Monomers ({len(cat.concatemer_3x_seq)} nt single-stranded concatemer)",
        f"Probe Binding Sites : {cat.repeating_probe_sites_count} binding sites detected across 3 repeats",
    ])
    
    if report.probe_dna:
        lines.extend([
            f"Detection Probe     : 5'-[FAM]-{report.probe_dna}-3' ({len(report.probe_dna)} nt)",
            f"Probe Tm / GC%      : {cat.probe_tm:.1f}°C  |  {cat.probe_gc:.1f}% (Binding dG: {cat.probe_dg:.2f} kcal/mol)",
            "-" * 90,
            "Concatemer Probe Hybridization Architecture (Repeated along concatemer):",
            "  Concatemer Strand   : 5' ... [ Repeat Unit 1 ] ─── [ Repeat Unit 2 ] ─── [ Repeat Unit 3 ] ... 3'",
            "                                     │                       │                       │",
            f"  Detection Probes    :        [FAM]-Probe             [FAM]-Probe             [FAM]-Probe",
            f"                              ({cat.probe_tm:.1f}°C)                 ({cat.probe_tm:.1f}°C)                 ({cat.probe_tm:.1f}°C)",
            "-" * 90,
        ])
        
    return "\n".join(lines)


def generate_qc_matrix_ascii(report: SimulationReport) -> str:
    """
    Renders cross-reactivity and biophysical QC matrix.
    """
    qc = report.qc_matrix
    tb = report.target_binding
    
    lines = [
        "=" * 90,
        "STAGE 5: IN SILICO BIOPHYSICAL QUALITY CONTROL & CROSS-REACTIVITY MATRIX",
        "=" * 90,
        f"Overall Assay Verdict : [{report.overall_verdict}]",
        "-" * 90,
        "Biophysical Checklist:",
        f"  1. Flush Ligation Nick (0 bp gap)        : [{'PASS' if tb.is_flush_nick else 'FAIL'}]",
        f"  2. Nick Junction -1 Base Complementarity : [{'PASS' if tb.nick_minus1_match else 'FAIL'}]",
        f"  3. Nick Junction +1 Base Complementarity : [{'PASS' if tb.nick_plus1_match else 'FAIL'}]",
        f"  4. Left Arm Tm >= 37°C ({tb.left_arm_tm:.1f}°C)        : [{'PASS' if tb.left_arm_tm >= report.reaction_temp else 'WARN'}]",
        f"  5. Right Arm Tm >= 37°C ({tb.right_arm_tm:.1f}°C)       : [{'PASS' if tb.right_arm_tm >= report.reaction_temp else 'WARN'}]",
        f"  6. Arm Tm Difference <= 5°C (ΔTm={tb.tm_diff:.1f}°C)   : [{'PASS' if tb.tm_diff <= 5.0 else 'WARN'}]",
        f"  7. RCA Primer Hairpin Check              : [{'PASS' if qc.primer_hairpin_found == 'No' else 'WARN'}] (Tm: {qc.primer_hairpin_tm:.1f}°C)",
        f"  8. Detection Probe Hairpin Check         : [{'PASS' if qc.probe_hairpin_found == 'No' else 'WARN'}] (Tm: {qc.probe_hairpin_tm:.1f}°C)",
        f"  9. Primer:Probe Cross-Dimer Check        : [{'PASS' if qc.primer_probe_dimer_dg > -3.0 else 'WARN'}] (dG: {qc.primer_probe_dimer_dg:.2f} kcal/mol)",
        "-" * 90,
    ]
    if qc.warnings:
        lines.append("Warnings Detected:")
        for w in qc.warnings:
            lines.append(f"  [!] {w}")
        lines.append("-" * 90)
    else:
        lines.append("Zero cross-dimer or thermodynamic interference detected. Probe and primer set fully orthogonal.")
        lines.append("-" * 90)
        
    return "\n".join(lines)


def render_full_terminal_report(report: SimulationReport) -> str:
    """Combines all stages into a complete terminal output string."""
    return "\n\n".join([
        generate_ascii_alignment(report),
        generate_circular_topology_ascii(report),
        generate_rca_and_detection_ascii(report),
        generate_qc_matrix_ascii(report),
    ])
