import glob
import os
import shutil
import subprocess
import tempfile

from pyPLP.thermodynamics.ct_parser import ct_to_dot_bracket
from pyPLP.thermodynamics.sir_graph import run_sir_graph_idt_plot

HAS_MFOLD = shutil.which("mfold") is not None
HAS_RNAPLOT = shutil.which("RNAplot") is not None


def eval_structure_mfold_detailed(
    seq_str: str,
    seq_id: str = "PLP",
    output_dir: str = None,
    temperature: float = 37.0,
    na_conc: float = 0.05,
    mg_conc: float = 0.0015,
    keep_ct_ps: bool = False,
):
    """Evaluates DNA minimum free energy (dG in kcal/mol) using mfold.

    Cleans up .ct/.ps based on keep_ct_ps.
    """
    if not HAS_MFOLD:
        return None, "N/A", "N/A"

    clean_seq = str(seq_str).strip().upper()
    seq_len = len(clean_seq)
    if seq_len < 4:
        return 0.0, "." * seq_len, "N/A"

    with tempfile.TemporaryDirectory() as tmp_dir:
        seq_filename = "query.seq"
        with open(os.path.join(tmp_dir, seq_filename), "w", encoding="ascii") as f:
            f.write(f">query\n{clean_seq}\n")

        cmd = [
            "mfold",
            f"SEQ={seq_filename}",
            "NA=DNA",
            "RUN_TYPE=text",
        ]

        subprocess.run(
            cmd,
            cwd=tmp_dir,
            capture_output=True,
            check=False,
            text=True,
            start_new_session=True,  # detaches from controlling tty (setsid)
        )

        # 1. Parse calculated dG safely with Latin-1 encoding
        dg = 0.0
        det_file = os.path.join(tmp_dir, "query.det")
        if os.path.exists(det_file):
            with open(det_file, "r", encoding="latin-1", errors="replace") as f:
                for line in f:
                    if "dG =" in line:
                        parts = line.replace("=", " ").split()
                        for i, token in enumerate(parts):
                            if token == "dG" and i + 1 < len(parts):
                                try:
                                    dg = float(parts[i + 1])
                                    break
                                except ValueError:
                                    pass
                        if dg != 0.0:
                            break

        if dg == 0.0:
            dg_file = os.path.join(tmp_dir, "query.dG")
            if os.path.exists(dg_file):
                with open(dg_file, "r", encoding="latin-1", errors="replace") as f:
                    lines = [l.strip() for l in f if l.strip()]
                    if lines:
                        try:
                            dg = float(lines[-1].split()[0])
                        except (ValueError, IndexError):
                            pass

        # 2. Locate generated .ct file
        ct_candidates = glob.glob(os.path.join(tmp_dir, "*.ct"))
        primary_ct = None
        for candidate in ct_candidates:
            if candidate.endswith("_1.ct"):
                primary_ct = candidate
                break
        if not primary_ct and ct_candidates:
            primary_ct = ct_candidates[0]

        dot_bracket = "." * seq_len
        has_pairs = False
        if primary_ct and os.path.exists(primary_ct):
            with open(primary_ct, "r", encoding="latin-1", errors="replace") as f:
                ct_lines = f.readlines()
                if len(ct_lines) > 1:
                    dot_bracket = ct_to_dot_bracket(ct_lines[1:], seq_len)
                    has_pairs = "(" in dot_bracket

        # 3. Export structure plots to output_dir
        plot_file_path = "N/A"
        if output_dir:
            abs_out_dir = os.path.abspath(output_dir)
            os.makedirs(abs_out_dir, exist_ok=True)
            dest_ct = os.path.join(abs_out_dir, f"{seq_id}.ct")

            if primary_ct and os.path.exists(primary_ct) and has_pairs:
                shutil.copyfile(primary_ct, dest_ct)
                plot_file_path = run_sir_graph_idt_plot(
                    dest_ct, abs_out_dir, seq_id, dg, keep_ct_ps=keep_ct_ps
                )
            else:
                with open(dest_ct, "w", encoding="ascii") as f_ct:
                    f_ct.write(f" {seq_len}  dG = {dg:.3f}  {seq_id}\n")
                    for idx, base in enumerate(clean_seq, start=1):
                        prev_i = idx - 1
                        next_i = idx + 1 if idx < seq_len else 0
                        f_ct.write(f" {idx:5d} {base}   {prev_i:5d} {next_i:5d}     0 {idx:5d}\n")

                if HAS_RNAPLOT:
                    svg_dest = os.path.join(abs_out_dir, f"{seq_id}_structure.svg")
                    plot_payload = f">{seq_id}\n{clean_seq}\n{dot_bracket}\n"
                    subprocess.run(
                        ["RNAplot", "--output-format=svg", "--layout-type=1"],
                        input=plot_payload,
                        text=True,
                        capture_output=True,
                        cwd=abs_out_dir,
                        check=False,
                    )
                    default_svg = os.path.join(abs_out_dir, f"{seq_id}_ss.svg")
                    if os.path.exists(default_svg):
                        shutil.move(default_svg, svg_dest)
                        plot_file_path = svg_dest

                if not keep_ct_ps and os.path.exists(dest_ct):
                    try:
                        os.remove(dest_ct)
                    except OSError:
                        pass

        return dg, dot_bracket, plot_file_path
