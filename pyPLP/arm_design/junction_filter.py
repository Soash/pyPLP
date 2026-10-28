"""pyPLP.arm_design.junction_filter: Quality filtering for padlock probe ligation junctions.

Evaluates terminal stacking energy at the nick and checks against nick self-dimerization
(target-independent probe circularization).
"""

from Bio.Seq import Seq


def evaluate_ligation_junction(
    left_arm: str,
    right_arm: str,
    min_gc_at_junction: int = 1,
    max_terminal_at_run: int = 2,
    check_nick_dimer: bool = True,
    max_nick_dimer_bp: int = 3,
) -> tuple[bool, str, dict]:
    """Evaluates the biochemical and structural quality of the padlock probe nick junction.

    The nick is formed by:
      - Position -1 (3'-OH side): right_arm[-1]
      - Position +1 (5'-PO4 side): left_arm[0]

    Args:
        left_arm: 5'->3' Left Arm sequence (5' end forms position +1 of the nick).
        right_arm: 5'->3' Right Arm sequence (3' end forms position -1 of the nick).
        min_gc_at_junction: Minimum G/C nucleotides required between right_arm[-1] and left_arm[0].
        max_terminal_at_run: Maximum permitted consecutive A/T bases at the 3' terminus of right_arm.
        check_nick_dimer: If True, checks for self-complementarity between the terminal ends.
        max_nick_dimer_bp: Maximum allowed consecutive Watson-Crick complementary bases between nick termini.

    Returns:
        (passes_filter, reason_str, metrics_dict)
    """
    l_arm = str(left_arm).strip().upper()
    r_arm = str(right_arm).strip().upper()

    if not l_arm or not r_arm:
        return False, "Empty arm sequence", {}

    # 1. Junction bases
    base_3prime = r_arm[-1]  # Position -1 (3'-OH)
    base_5prime = l_arm[0]   # Position +1 (5'-PO4)

    gc_count = sum(1 for b in (base_3prime, base_5prime) if b in ("G", "C"))

    # 2. Terminal A/T run at 3' end of right_arm
    at_run = 0
    for b in reversed(r_arm):
        if b in ("A", "T"):
            at_run += 1
        else:
            break

    # 3. Nick self-dimer check (terminal 5 nt of each end)
    nick_dimer_found = False
    max_matched_bp = 0
    if check_nick_dimer:
        # 3' end of right_arm (5 nt) vs 5' end of left_arm (5 nt)
        term_r = r_arm[-5:] if len(r_arm) >= 5 else r_arm
        term_l = l_arm[:5] if len(l_arm) >= 5 else l_arm

        # If term_r is complementary to term_l in antiparallel orientation:
        # term_l 5'->3' vs term_r 3'->5'
        term_r_rev = term_r[::-1]  # 3' to 5'
        comp_map = {"A": "T", "T": "A", "C": "G", "G": "C"}
        
        # Check contiguous complementary matches between the two ends
        consec = 0
        for b_l, b_r in zip(term_l, term_r_rev):
            if comp_map.get(b_l) == b_r:
                consec += 1
                if consec > max_matched_bp:
                    max_matched_bp = consec
            else:
                consec = 0

        if max_matched_bp > max_nick_dimer_bp:
            nick_dimer_found = True

    metrics = {
        "Junction_3prime_Base": base_3prime,
        "Junction_5prime_Base": base_5prime,
        "Junction_GC": gc_count,
        "Terminal_3prime_AT_Run": at_run,
        "Nick_Dimer_Max_BP": max_matched_bp,
    }

    if gc_count < min_gc_at_junction:
        return False, f"Low junction GC ({gc_count} < {min_gc_at_junction})", metrics

    if at_run > max_terminal_at_run:
        return False, f"Long terminal A/T run ({at_run} > {max_terminal_at_run})", metrics

    if nick_dimer_found:
        return False, f"Nick self-dimerization detected ({max_matched_bp} bp)", metrics

    return True, "Pass", metrics
