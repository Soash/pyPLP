from collections import Counter


def get_column_consensus_and_conservation(column: str, ignore_gaps: bool = False):
    """Returns the most frequent nucleotide (A/C/G/T) and its conservation percentage in an alignment column."""
    col = column.upper()
    total_len = len(col)
    if total_len == 0:
        return "-", 0.0

    nuc_bases = [b for b in col if b in "ACGT"]
    if not nuc_bases:
        return "-", 0.0

    counts = Counter(nuc_bases)
    most_common_base, most_common_count = counts.most_common(1)[0]

    denom = len(col.replace("-", "").replace(".", "")) if ignore_gaps else total_len
    if denom == 0:
        return "-", 0.0

    score = (most_common_count / denom) * 100.0
    return most_common_base, score
