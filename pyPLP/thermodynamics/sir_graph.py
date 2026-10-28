import glob
import os
import shutil
import subprocess

from pyPLP.thermodynamics.ps_converter import convert_ps_to_png


def run_sir_graph_idt_plot(ct_path: str, output_dir: str, seq_id: str, dG: float, keep_ct_ps: bool = False) -> str:
    """Executes sir_graph without -c to preserve 5' and 3' terminal open markers."""
    sir_bin = shutil.which("sir_graph") or os.path.expanduser("~/miniconda3/bin/sir_graph")
    if not (sir_bin and os.path.exists(ct_path)):
        return "N/A"

    abs_out_dir = os.path.abspath(output_dir)
    ct_name = os.path.basename(ct_path)
    title_str = f"dG = {dG:.3f} kcal/mol  {seq_id}"

    cmd = [
        sir_bin,
        "-p",
        "-s",
        "-D",
        "-lab", "10",
        "-t", title_str,
        "-o", seq_id,
        ct_name,
    ]

    subprocess.run(
        cmd,
        cwd=abs_out_dir,
        capture_output=True,
        check=False,
    )

    ps_file = os.path.join(abs_out_dir, f"{seq_id}.ps")
    if not os.path.exists(ps_file):
        candidates = glob.glob(os.path.join(abs_out_dir, f"{seq_id}*.ps"))
        if candidates:
            ps_file = candidates[0]

    final_plot = "N/A"
    if os.path.exists(ps_file):
        target_png = os.path.join(abs_out_dir, f"{seq_id}_structure.png")
        result = convert_ps_to_png(ps_file, target_png)
        if result.endswith(".png") and os.path.exists(result):
            final_plot = result
            if not keep_ct_ps:
                try:
                    os.remove(ps_file)
                except OSError:
                    pass
        else:
            final_plot = ps_file

    if not keep_ct_ps and os.path.exists(ct_path):
        try:
            os.remove(ct_path)
        except OSError:
            pass

    return final_plot
