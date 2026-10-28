import os
import shutil
import subprocess

HAS_RNAPLOT = shutil.which("RNAplot") is not None


def render_vienna_svg(
    seq_id: str,
    full_seq: str,
    dot_bracket: str,
    vienna_dir: str,
) -> str:
    """Renders a 2D secondary structure diagram using ViennaRNA RNAplot in SVG format."""
    svg_file = os.path.join(vienna_dir, f"{seq_id}_structure.svg")
    if HAS_RNAPLOT and dot_bracket not in ("N/A", ""):
        plot_payload = f">{seq_id}\n{full_seq}\n{dot_bracket}\n"
        subprocess.run(
            ["RNAplot", "--output-format=svg", "--layout-type=1"],
            input=plot_payload,
            text=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            cwd=vienna_dir,
            check=False,
        )
        rna_svg = os.path.join(vienna_dir, f"{seq_id}_ss.svg")
        if os.path.exists(rna_svg):
            shutil.move(rna_svg, svg_file)
            return svg_file
    return "N/A"
