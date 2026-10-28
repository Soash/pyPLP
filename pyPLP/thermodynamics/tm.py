import primer3


def calculate_tm(seq_str: str, **kwargs) -> float:
    """Calculates nearest-neighbor melting temperature (Tm) using Primer3 defaults:

    mv_conc=50.0 mM, dv_conc=1.5 mM, dntp_conc=0.6 mM, dna_conc=50.0 nM, tm_method='santalucia'.
    Optional keyword arguments can override any parameter.
    """
    clean_seq = str(seq_str).strip().upper()
    if not clean_seq:
        return 0.0

    tm = primer3.calc_tm(clean_seq, **kwargs)
    return round(float(tm), 2)


# Alias for backward compatibility
calculate_tm_oligoanalyzer = calculate_tm
