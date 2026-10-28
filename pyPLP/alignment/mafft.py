from pathlib import Path
import subprocess
import sys


def run_mafft(
    input_fasta: Path | str,
    output_fasta: Path | str,
    mafft_env: str = "mafft_env",
    mafft_mode: str = "fast",
    mafft_threads: int | str = 8,
) -> bool:
    """Executes MAFFT alignment via conda environment."""
    in_path = Path(input_fasta).resolve()
    out_path = Path(output_fasta).resolve()
    mode = str(mafft_mode).lower()
    threads = str(mafft_threads)

    mafft_args = (
        ["mafft", "--thread", threads, "--localpair", "--maxiterate", "1000", str(in_path)]
        if mode == "slow"
        else ["mafft", "--thread", threads, "--auto", str(in_path)]
    )
    full_cmd = ["conda", "run", "-n", mafft_env, "--no-capture-output"] + mafft_args
    print(f"[MAFFT] Executing MAFFT in conda env '{mafft_env}' ({mode} mode, {threads} threads)...")

    with open(out_path, "w") as out_aln:
        res = subprocess.run(full_cmd, stdout=out_aln, text=True)
        if res.returncode != 0:
            print(f"[Error] MAFFT run failed in conda environment '{mafft_env}'.", file=sys.stderr)
            return False
    return True
