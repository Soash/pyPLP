"""pyPLP.multiplex.cross_dimer: Pairwise thermodynamic cross-dimerization calculations.
"""

import primer3


def calc_heterodimer_dg(
    seq1: str,
    seq2: str,
    temp_c: float = 37.0,
    min_dimer_tm: float = 15.0,
    mv_conc: float = 50.0,
    dv_conc: float = 1.5,
    dntp_conc: float = 0.6,
    dna_conc: float = 50.0,
) -> float:
    """Calculates heterodimer minimum free energy (dG in kcal/mol at temp_c)
    between two DNA oligonucleotides using Primer3 SantaLucia nearest-neighbor model.
    Duplexes with Tm below min_dimer_tm are unstable at assay temperature and return 0.0.
    """
    s1 = str(seq1).strip().upper()
    s2 = str(seq2).strip().upper()
    if not s1 or not s2:
        return 0.0

    try:
        res = primer3.calc_heterodimer(
            s1,
            s2,
            temp_c=temp_c,
            mv_conc=mv_conc,
            dv_conc=dv_conc,
            dntp_conc=dntp_conc,
            dna_conc=dna_conc,
        )
        if res.structure_found and res.tm >= min_dimer_tm and res.dg != 0.0:
            return round(float(res.dg) / 1000.0, 2)
        return 0.0
    except Exception:
        return 0.0


def calc_homodimer_dg(
    seq: str,
    temp_c: float = 37.0,
    min_dimer_tm: float = 15.0,
    mv_conc: float = 50.0,
    dv_conc: float = 1.5,
    dntp_conc: float = 0.6,
    dna_conc: float = 50.0,
) -> float:
    """Calculates homodimer self-folding free energy (dG in kcal/mol at temp_c)
    for a single DNA oligonucleotide.
    Duplexes with Tm below min_dimer_tm are unstable at assay temperature and return 0.0.
    """
    s = str(seq).strip().upper()
    if not s:
        return 0.0

    try:
        res = primer3.calc_homodimer(
            s,
            temp_c=temp_c,
            mv_conc=mv_conc,
            dv_conc=dv_conc,
            dntp_conc=dntp_conc,
            dna_conc=dna_conc,
        )
        if res.structure_found and res.tm >= min_dimer_tm and res.dg != 0.0:
            return round(float(res.dg) / 1000.0, 2)
        return 0.0
    except Exception:
        return 0.0


def eval_oligo_orthogonality(
    named_oligos: dict[str, str],
    max_cross_dg: float = -2.0,
    temp_c: float = 37.0,
) -> tuple[bool, list[dict]]:
    """Evaluates all pairwise interactions across a pool of named oligonucleotides.

    Flags any interaction where dG < max_cross_dg (e.g. -3.5 kcal/mol < -2.0 kcal/mol).
    Returns (is_orthogonal, violations_list).
    """
    names = list(named_oligos.keys())
    violations = []

    for i in range(len(names)):
        n1 = names[i]
        s1 = named_oligos[n1]

        # Check homodimer
        self_dg = calc_homodimer_dg(s1, temp_c=temp_c)
        if self_dg < max_cross_dg:
            violations.append({
                "oligo1": n1,
                "oligo2": n1,
                "type": "homodimer",
                "dg": self_dg,
                "cutoff": max_cross_dg,
            })

        for j in range(i + 1, len(names)):
            n2 = names[j]
            s2 = named_oligos[n2]
            cross_dg = calc_heterodimer_dg(s1, s2, temp_c=temp_c)
            if cross_dg < max_cross_dg:
                violations.append({
                    "oligo1": n1,
                    "oligo2": n2,
                    "type": "heterodimer",
                    "dg": cross_dg,
                    "cutoff": max_cross_dg,
                })

    is_orthogonal = (len(violations) == 0)
    return is_orthogonal, violations
