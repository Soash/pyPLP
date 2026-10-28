import os
import shutil
import subprocess


def convert_ps_to_png(ps_file: str, output_png: str) -> str:
    """Converts PostScript to centered PNG via Ghostscript with -dEPSCrop."""
    gs_bin = shutil.which("gs")
    if gs_bin and os.path.exists(ps_file):
        cmd = [
            gs_bin,
            "-dSAFER",
            "-dBATCH",
            "-dNOPAUSE",
            "-dEPSCrop",
            "-sDEVICE=png16m",
            "-r300",
            f"-sOutputFile={output_png}",
            ps_file,
        ]
        env = os.environ.copy()
        env["TMPDIR"] = "/tmp"

        res = subprocess.run(
            cmd,
            env=env,
            capture_output=True,
            check=False,
        )
        if res.returncode == 0 and os.path.exists(output_png):
            return output_png

    return ps_file
