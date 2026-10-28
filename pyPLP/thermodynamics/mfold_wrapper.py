from pyPLP.thermodynamics.mfold_detailed import eval_structure_mfold_detailed


def eval_structure_mfold(
    seq_str: str,
    temperature: float = 37.0,
    na_conc: float = 0.05,
    mg_conc: float = 0.0015,
) -> float | None:
    """Backward-compatible wrapper returning minimum free energy (dG in kcal/mol)."""
    dg, _, _ = eval_structure_mfold_detailed(
        seq_str,
        temperature=temperature,
        na_conc=na_conc,
        mg_conc=mg_conc,
    )
    return dg
