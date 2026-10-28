import os
import shutil
import subprocess

from pyPLP.arm_design.dimer_eval import get_vienna_dna_parameter_path
from pyPLP.thermodynamics.mfold_detailed import eval_structure_mfold_detailed

HAS_RNAFOLD = shutil.which("RNAfold") is not None
HAS_MFOLD = shutil.which("mfold") is not None


def eval_plp_structure(
    seq_str: str,
    folding_engine: str = "RNAfold",
    temperature: float = 37.0,
) -> tuple[str, float]:
    """Evaluates full PLP secondary structure using RNAfold (with DNA parameters) or mfold.

    Returns (dot_bracket, mfe_dg).
    """
    clean_seq = str(seq_str).strip().upper()

    if folding_engine.lower() == "mfold" and HAS_MFOLD:
        dg, dot_bracket, _ = eval_structure_mfold_detailed(
            clean_seq,
            seq_id="screen",
            output_dir=None,
            temperature=temperature,
        )
        return dot_bracket, (dg if dg is not None else 0.0)

    # Default to RNAfold
    if not HAS_RNAFOLD:
        return "N/A", 0.0

    cmd = ["RNAfold", "--noconv", "-T", str(temperature), "--noPS"]
    dna_param_path = get_vienna_dna_parameter_path()
    if dna_param_path and os.path.exists(dna_param_path):
        cmd.extend(["-P", dna_param_path])

    res = subprocess.run(
        cmd,
        input=f">q\n{clean_seq}\n",
        text=True,
        capture_output=True,
        check=False,
    )
    if res.returncode != 0:
        return "N/A", 0.0

    lines = res.stdout.strip().split("\n")
    if len(lines) >= 3:
        parts = lines[2].split()
        dot_bracket = parts[0]
        val_str = lines[2].split("(")[-1].replace(")", "").replace("=", "").strip()
        try:
            mfe_dg = float(val_str)
        except ValueError:
            mfe_dg = 0.0
        return dot_bracket, mfe_dg

    return "N/A", 0.0
