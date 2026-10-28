import re


def get_max_homopolymer(seq_str: str):
    """Finds the longest run of identical consecutive bases in a sequence."""
    clean_seq = str(seq_str).strip().upper()
    if not clean_seq:
        return 0, ""

    matches = re.findall(r"((.)\2*)", clean_seq)
    if not matches:
        return 0, ""

    longest = max(matches, key=lambda x: len(x[0]))
    return len(longest[0]), longest[0]
