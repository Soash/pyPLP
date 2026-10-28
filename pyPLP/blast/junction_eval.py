"""pyPLP.blast.junction_eval: Ligase footprint and 3' nick mismatch evaluation for BLAST hits.

Evaluates whether off-target BLAST hits align across the ligation junction and whether
the critical 3' terminal nucleotide (position -1 from the nick) is matched or mismatched.
"""


def evaluate_hit_ligation_competence(
    qstart: int,
    qend: int,
    qseq: str,
    sseq: str,
    arm_len: int,
) -> tuple[bool, str]:
    """Determines whether a BLAST off-target alignment is enzymatically competent for ligation.

    In a contiguous candidate footprint (2 * arm_len):
      - 1 .. arm_len: Left Target (binds Left Arm, 5'-phosphate at coordinate arm_len)
      - arm_len + 1 .. 2*arm_len: Right Target (binds Right Arm, 3'-OH at coordinate arm_len + 1)
      - Ligation nick: between coordinate arm_len and arm_len + 1.
      - Position -1 (critical 3' base): coordinate arm_len + 1.

    Args:
        qstart: 1-based start position of the alignment in the query footprint.
        qend: 1-based end position of the alignment in the query footprint.
        qseq: Aligned query sequence (may contain '-' for deletions).
        sseq: Aligned subject sequence (may contain '-' for insertions).
        arm_len: Length of one arm (half of the full footprint).

    Returns:
        (is_ligation_competent, classification_reason)
    """
    nick_3prime_pos = arm_len + 1
    nick_5prime_pos = arm_len

    # If alignment does not span the critical 3' junction base
    if qstart > nick_3prime_pos or qend < nick_3prime_pos:
        return False, "Alignment does not span 3' nick position"

    # Walk through alignment string to find the column for nick_3prime_pos
    curr_qpos = qstart
    base_3p_match = False
    base_5p_match = False

    for q_char, s_char in zip(qseq.upper(), sseq.upper()):
        if q_char != "-":
            if curr_qpos == nick_3prime_pos:
                # 3' terminus of Right Arm
                if s_char != "-" and q_char == s_char:
                    base_3p_match = True
                else:
                    return False, f"Mismatch at 3' nick position ({q_char} vs {s_char})"

            if curr_qpos == nick_5prime_pos:
                # 5' terminus of Left Arm
                if s_char != "-" and q_char == s_char:
                    base_5p_match = True

            curr_qpos += 1

    if base_3p_match:
        if base_5p_match:
            return True, "Intact nick junction (matched at both 3' and 5' nick bases)"
        return True, "Matched at critical 3' nick base"

    return False, "Unmatched 3' nick base"
