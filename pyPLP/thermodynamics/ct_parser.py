def ct_to_dot_bracket(ct_lines, seq_length: int) -> str:
    """Converts mfold Connectivity Table (.ct) format lines to dot-bracket notation."""
    pairs = [0] * (seq_length + 1)
    for line in ct_lines:
        tokens = line.strip().split()
        if len(tokens) >= 6 and tokens[0].isdigit() and tokens[4].isdigit():
            i = int(tokens[0])
            j = int(tokens[4])
            if i <= seq_length:
                pairs[i] = j

    db = []
    for i in range(1, seq_length + 1):
        j = pairs[i]
        if j == 0:
            db.append(".")
        elif i < j:
            db.append("(")
        else:
            db.append(")")
    return "".join(db)
