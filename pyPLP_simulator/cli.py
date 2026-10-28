#!/usr/bin/env python3
"""
pyPLP_simulator.cli: Command-line interface for the in silico padlock probe reaction simulator.
"""

import os
import sys
import json
import argparse
from pathlib import Path
from typing import Optional

from pyPLP_simulator.core import SimulationReport, clean_sequence
from pyPLP_simulator.visualizer import render_full_terminal_report
from pyPLP_simulator.html_report import generate_interactive_html


def load_sequence_or_string(val: str, name: str) -> str:
    """Loads sequence from a file path if it exists, or returns raw string."""
    if not val:
        return ""
    stripped = str(val).strip()
    if os.path.isfile(stripped):
        try:
            with open(stripped, "r", encoding="utf-8") as f:
                content = f.read()
                # Check FASTA
                lines = [l.strip() for l in content.splitlines() if l.strip()]
                seq_lines = [l for l in lines if not l.startswith(">")]
                raw_seq = "".join(seq_lines)
                return clean_sequence(raw_seq, name)
        except Exception as e:
            raise ValueError(f"Could not read sequence from file '{stripped}': {e}")
    return clean_sequence(stripped, name)


def run_simulation(
    target: str,
    plp: str,
    rca_primer: Optional[str] = None,
    probe: Optional[str] = None,
    output_html: Optional[str] = None,
    output_json: Optional[str] = None,
    job_name: Optional[str] = None,
    reaction_temp: float = 37.0,
    is_rna: bool = False,
    quiet: bool = False,
) -> SimulationReport:
    """Executes the full in silico molecular reaction simulation and exports reports."""
    clean_target = load_sequence_or_string(target, "Target")
    clean_plp = load_sequence_or_string(plp, "Padlock Probe (PLP)")
    clean_primer = load_sequence_or_string(rca_primer, "RCA Primer") if rca_primer else None
    clean_probe = load_sequence_or_string(probe, "Detection Probe") if probe else None

    # Determine destination directory
    job_dir = "."
    if job_name:
        raw_job = str(job_name).strip()
        job_path = Path(raw_job)
        if not job_path.is_absolute():
            parts = job_path.parts
            if not parts or parts[0] != "results":
                job_path = Path("results") / job_path
        job_dir = str(job_path.resolve())
        os.makedirs(job_dir, exist_ok=True)

    # Run core simulation
    report = SimulationReport(
        target_seq=clean_target,
        plp_seq=clean_plp,
        rca_primer_seq=clean_primer,
        probe_seq=clean_probe,
        reaction_temp=reaction_temp,
        is_rna=is_rna,
    )

    # Print terminal diagrams
    if not quiet:
        term_output = render_full_terminal_report(report)
        print(term_output)

    # Export interactive HTML visualizer
    final_html_path = None
    if output_html:
        if os.path.isabs(output_html) or os.path.dirname(output_html):
            final_html_path = os.path.abspath(output_html)
        else:
            final_html_path = os.path.join(job_dir, output_html)
    elif job_name:
        final_html_path = os.path.join(job_dir, "simulation_report.html")
    else:
        final_html_path = os.path.abspath("simulation_report.html")

    if final_html_path:
        os.makedirs(os.path.dirname(final_html_path), exist_ok=True)
        generate_interactive_html(report, final_html_path)
        if not quiet:
            print(f"\n[HTML Visualizer] Interactive molecular report saved to:")
            print(f"  -> file://{final_html_path}\n")

    # Export machine-readable JSON if requested
    if output_json:
        if os.path.isabs(output_json) or os.path.dirname(output_json):
            final_json_path = os.path.abspath(output_json)
        else:
            final_json_path = os.path.join(job_dir, output_json)
        os.makedirs(os.path.dirname(final_json_path), exist_ok=True)
        
        json_data = {
            "verdict": report.overall_verdict,
            "target_type": "RNA" if report.is_rna else "DNA",
            "reaction_temp_c": report.reaction_temp,
            "left_arm": {
                "sequence": report.target_binding.left_arm_seq,
                "length": report.target_binding.left_arm_len,
                "tm": report.target_binding.left_arm_tm,
                "dg": report.target_binding.left_arm_dg,
            },
            "right_arm": {
                "sequence": report.target_binding.right_arm_seq,
                "length": report.target_binding.right_arm_len,
                "tm": report.target_binding.right_arm_tm,
                "dg": report.target_binding.right_arm_dg,
            },
            "ligation_nick": {
                "is_flush": report.target_binding.is_flush_nick,
                "gap_bp": report.target_binding.gap_bp,
                "minus_1_match": report.target_binding.nick_minus1_match,
                "plus_1_match": report.target_binding.nick_plus1_match,
            },
            "rca_priming": {
                "primer_tm": report.rca_priming.primer_tm if report.rca_priming else None,
                "primer_dg": report.rca_priming.primer_dg if report.rca_priming else None,
            },
            "concatemer": {
                "monomer_length": report.concatemer.unit_repeat_len,
                "probe_sites_in_3_repeats": report.concatemer.repeating_probe_sites_count,
                "probe_tm": report.concatemer.probe_tm,
            },
            "qc_cross_reactivity": {
                "primer_probe_dimer_dg": report.qc_matrix.primer_probe_dimer_dg,
                "primer_self_dimer_dg": report.qc_matrix.primer_self_dimer_dg,
                "probe_self_dimer_dg": report.qc_matrix.probe_self_dimer_dg,
                "warnings": report.qc_matrix.warnings,
            },
        }
        with open(final_json_path, "w", encoding="utf-8") as jf:
            json.dump(json_data, jf, indent=2)

    return report


def main():
    parser = argparse.ArgumentParser(
        description="pyPLP_simulator: In silico molecular reaction simulator and visualizer for Padlock Probes (PLP), RCA Primers, and Detection Probes."
    )
    parser.add_argument("-t", "--target", help="Target sequence (5' -> 3', RNA or DNA) or path to FASTA file")
    parser.add_argument("-p", "--plp", help="Padlock Probe (PLP) sequence (5' -> 3') or path to sequence file")
    parser.add_argument("-r", "--rca-primer", default=None, help="RCA Primer sequence (5' -> 3')")
    parser.add_argument("-d", "--probe", default=None, help="Detection Probe sequence (5' -> 3', e.g. 5'-FAM labeled)")
    parser.add_argument("-T", "--temp", type=float, default=37.0, help="Isothermal reaction temperature in °C (default: 37.0)")
    parser.add_argument("--rna", action="store_true", help="Explicitly specify that target is RNA")
    parser.add_argument("-j", "--job-name", default=None, help="Output job directory name under results/")
    parser.add_argument("-o", "--output-html", default=None, help="Output HTML report filename")
    parser.add_argument("--json", default=None, help="Optional JSON output filepath")
    parser.add_argument("-q", "--quiet", action="store_true", help="Suppress terminal printout and only write HTML/JSON files")
    parser.add_argument("-i", "--interactive", action="store_true", help="Launch interactive terminal prompt to enter sequences")
    parser.add_argument("-w", "--web", action="store_true", help="Launch interactive animated web app simulator in browser")
    parser.add_argument("--port", type=int, default=None, help="Custom port for local web server (default: auto 8080)")

    args = parser.parse_args()

    # Launch Web Simulator if --web is requested
    if args.web:
        from pyPLP_simulator.web_server import run_web_simulator
        run_web_simulator(port=args.port)
        sys.exit(0)

    # Interactive mode if requested or if mandatory arguments are missing
    if args.interactive or (not args.target and not args.plp):
        print("=" * 80)
        print("pyPLP_simulator: Interactive In Silico Molecular Assay Visualizer")
        print("=" * 80)
        try:
            target_in = input("Enter Target Sequence (DNA or RNA 5'->3'): ").strip()
            plp_in = input("Enter Padlock Probe (PLP) Sequence (5'->3'): ").strip()
            primer_in = input("Enter RCA Primer Sequence (optional, press Enter to skip): ").strip() or None
            probe_in = input("Enter Detection Probe Sequence (optional, press Enter to skip): ").strip() or None
            temp_str = input("Reaction Temperature in °C (default 37.0): ").strip()
            temp_in = float(temp_str) if temp_str else 37.0
            job_in = input("Job / Output Directory Name (optional): ").strip() or "plp_simulation"
            
            run_simulation(
                target=target_in,
                plp=plp_in,
                rca_primer=primer_in,
                probe=probe_in,
                job_name=job_in,
                reaction_temp=temp_in,
                output_html="simulation_report.html",
            )
            sys.exit(0)
        except (KeyboardInterrupt, EOFError):
            print("\nAborted.")
            sys.exit(1)

    if not args.target or not args.plp:
        parser.error("Both --target (-t) and --plp (-p) sequences are required (or use --interactive).")

    try:
        run_simulation(
            target=args.target,
            plp=args.plp,
            rca_primer=args.rca_primer,
            probe=args.probe,
            output_html=args.output_html,
            output_json=args.json,
            job_name=args.job_name,
            reaction_temp=args.temp,
            is_rna=args.rna,
            quiet=args.quiet,
        )
    except Exception as err:
        print(f"\n[Simulation Error] {err}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
