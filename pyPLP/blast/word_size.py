def get_optimal_word_size(arm_len: int) -> int:
    """Dynamically scales word_size based on arm length to maximize sensitivity

    while optimizing memory and runtime across CPU cores.
    """
    if arm_len < 15:
        return 7  # High sensitivity for 12-14 nt short queries
    elif arm_len < 18:
        return 8  # Balanced for 15-17 nt
    else:
        return 10  # Optimal speed/sensitivity for 18-20 nt probes
