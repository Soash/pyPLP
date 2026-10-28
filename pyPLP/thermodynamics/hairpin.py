import primer3


def check_hairpin(seq_str: str, **kwargs):
    """Calculates unimolecular hairpin stability using Primer3 nearest-neighbor engine.

    Primer3's internal thal C library enforces a maximum sequence length of 60 bp.
    For sequences <= 60 bp, calc_hairpin is evaluated directly.
    For sequences > 60 bp (e.g. long backbones), a 60 bp sliding window scans the
    molecule to detect any localized stem-loop formations, reporting the most stable
    structure (highest Tm, lowest dG).
    """
    clean_seq = str(seq_str).strip().upper()
    seq_len = len(clean_seq)
    if seq_len < 4:
        return "No", 0.0, 0.0

    if seq_len <= 60:
        try:
            hp = primer3.calc_hairpin(clean_seq, **kwargs)
            hairpin_found = "Yes" if hp.structure_found else "No"
            hairpin_tm = round(float(hp.tm), 2) if hp.structure_found else 0.0
            hairpin_dg = round(float(hp.dg) / 1000.0, 2) if hp.structure_found else 0.0
            return hairpin_found, hairpin_tm, hairpin_dg
        except Exception:
            pass

    # Sliding 60 bp window for sequences > 60 bp
    window_size = 60
    max_tm = 0.0
    min_dg = 0.0
    found = False

    for i in range(seq_len - window_size + 1):
        sub_seq = clean_seq[i : i + window_size]
        try:
            hp = primer3.calc_hairpin(sub_seq, **kwargs)
            if hp.structure_found:
                found = True
                if hp.tm > max_tm:
                    max_tm = float(hp.tm)
                dg_val = float(hp.dg) / 1000.0
                if dg_val < min_dg:
                    min_dg = dg_val
        except Exception:
            continue

    return ("Yes" if found else "No"), round(max_tm, 2), round(min_dg, 2)
