import random


def parse_base_pairs(dot_bracket: str) -> list[tuple[int, int]]:
    """Extracts paired base index tuples (i, j) from dot-bracket notation."""
    stack = []
    pairs = []
    for idx, ch in enumerate(dot_bracket):
        if ch == "(":
            stack.append(idx)
        elif ch == ")":
            if stack:
                pairs.append((stack.pop(), idx))
    return pairs


def break_stem_mutations_aggressive(
    seq_str: str,
    dot_bracket: str,
    len_l: int,
    len_r: int,
    max_homo: int = 4,
    protected_indices: set[int] | None = None,
) -> tuple[str, bool]:
    """Identifies base pairs involving the backbone region and mutates backbone
    nucleotides to non-pairing bases while respecting homopolymer limits and protected regions.
    """
    pairs = parse_base_pairs(dot_bracket)
    total_len = len(seq_str)
    bb_start = len_l
    bb_end = total_len - len_r
    seq_list = list(seq_str)

    if not pairs:
        return seq_str, False

    pairing_rules = {
        "A": ["T", "U"],
        "T": ["A", "G"],
        "U": ["A", "G"],
        "C": ["G"],
        "G": ["C", "T", "U"],
    }

    mutated_any = False
    random.shuffle(pairs)

    for i, j in pairs:
        mutable_choices = []
        if bb_start <= i < bb_end and (protected_indices is None or i not in protected_indices):
            mutable_choices.append((i, j))
        if bb_start <= j < bb_end and (protected_indices is None or j not in protected_indices):
            mutable_choices.append((j, i))

        if not mutable_choices:
            continue

        target_idx, partner_idx = random.choice(mutable_choices)
        orig_base = seq_list[target_idx]
        partner_base = seq_list[partner_idx]

        forbidden = set(pairing_rules.get(partner_base, []) + [orig_base])
        candidates = [b for b in ["A", "C", "G", "T"] if b not in forbidden]
        random.shuffle(candidates)

        for cand in candidates:
            temp_seq = seq_list[:]
            temp_seq[target_idx] = cand
            window = "".join(temp_seq[max(0, target_idx - 3) : target_idx + 4])
            if cand * max_homo not in window:
                seq_list[target_idx] = cand
                mutated_any = True
                break

    return "".join(seq_list), mutated_any
