import os
import shutil
import subprocess
import sys

from pyPLP.thermodynamics.mfold_wrapper import eval_structure_mfold

HAS_RNAFOLD = shutil.which("RNAfold") is not None


def get_vienna_dna_parameter_path() -> str | None:
    """Finds the path to dna_mathews2004.par in ViennaRNA installation."""
    candidates = [
        os.path.join(sys.prefix, "share", "ViennaRNA", "dna_mathews2004.par"),
        os.path.expanduser("~/miniconda3/share/ViennaRNA/dna_mathews2004.par"),
        "/usr/share/ViennaRNA/dna_mathews2004.par",
        "/usr/local/share/ViennaRNA/dna_mathews2004.par",
    ]
    for path in candidates:
        if os.path.isfile(path):
            return path
    return "dna_mathews2004.par"


def eval_arm_dimer_vienna(left_arm: str, right_arm: str, temperature: float = 37.0) -> tuple[str, float]:
    """Evaluates cross-hybridization (heterodimer) between the two probe arms

    using ViennaRNA (RNAfold) with native SantaLucia/Mathews DNA parameters (dna_mathews2004.par).
    Returns (vienna_db, vienna_dg).
    """
    seq_str = left_arm + right_arm
    vienna_db = "N/A"
    vienna_dg = 0.0

    if HAS_RNAFOLD:
        cmd = ["RNAfold"]
        dna_par = get_vienna_dna_parameter_path()
        if dna_par:
            cmd.extend(["-P", dna_par])
        cmd.extend(["--noconv", "-T", str(temperature), "--noPS"])

        res = subprocess.run(
            cmd,
            input=f">dimer\n{seq_str}\n",
            text=True,
            capture_output=True,
            check=False,
        )
        if res.returncode == 0:
            lines = res.stdout.strip().split("\n")
            if len(lines) >= 3:
                vienna_db = lines[2].split()[0]
                val_str = (
                    lines[2].split("(")[-1].replace(")", "").replace("=", "").strip()
                )
                try:
                    vienna_dg = float(val_str)
                except ValueError:
                    vienna_dg = 0.0

    return vienna_db, vienna_dg


def eval_arm_dimer_mfold(left_arm: str, right_arm: str, temperature: float = 37.0) -> float | None:
    """Evaluates cross-hybridization (heterodimer) between the two probe arms

    using mfold (SantaLucia DNA parameters).
    Returns mfold_dg.
    """
    seq_str = left_arm + right_arm
    return eval_structure_mfold(seq_str, temperature=temperature)


def eval_arm_dimer_dual(left_arm: str, right_arm: str, temperature: float = 37.0):
    """Evaluates cross-hybridization (heterodimer) between the two probe arms

    using both ViennaRNA (with DNA parameters) and mfold (DNA model).
    Returns (vienna_db, vienna_dg, mfold_dg).
    """
    vienna_db, vienna_dg = eval_arm_dimer_vienna(left_arm, right_arm, temperature=temperature)
    mfold_dg = eval_arm_dimer_mfold(left_arm, right_arm, temperature=temperature)
    return vienna_db, vienna_dg, mfold_dg
