def find_contiguous_blocks(conserved_indices, min_block_len_bp):
    """Groups consecutive 0-based indices into contiguous (start, end) intervals

    that meet or exceed min_block_len_bp.
    """
    if not conserved_indices:
        return []

    blocks = []
    block_start = conserved_indices[0]
    prev = conserved_indices[0]

    for idx in conserved_indices[1:]:
        if idx == prev + 1:
            prev = idx
        else:
            if (prev - block_start + 1) >= min_block_len_bp:
                blocks.append((block_start, prev))
            block_start = idx
            prev = idx

    if (prev - block_start + 1) >= min_block_len_bp:
        blocks.append((block_start, prev))

    return blocks
