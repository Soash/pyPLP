from pyPLP.thermodynamics.tm import calculate_tm


def compute_metrics(seq_obj):
    """Computes GC content and standard Primer3 default melting temperature (Tm)."""
    seq_str = str(seq_obj).strip().upper()
    if not seq_str:
        return 0.0, 0.0

    gc_count = seq_str.count("G") + seq_str.count("C")
    gc_percent = round((gc_count / len(seq_str)) * 100.0, 2)
    tm = calculate_tm(seq_str)

    return gc_percent, tm
