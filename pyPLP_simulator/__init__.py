"""
pyPLP_simulator: In Silico Padlock Probe (PLP) & RCA Reaction Simulator & Visualizer.
"""

from pyPLP_simulator.core import (
    SimulationReport,
    TargetBindingResult,
    BackboneAnnotationResult,
    RCAPrimingResult,
    ConcatemerResult,
    QCMatrixResult,
    find_best_arm_mapping,
    annotate_backbone_regions,
    simulate_rca_priming,
    generate_rca_concatemer,
    run_cross_reactivity_qc,
)
from pyPLP_simulator.visualizer import render_full_terminal_report
from pyPLP_simulator.html_report import generate_interactive_html
from pyPLP_simulator.cli import run_simulation
from pyPLP_simulator.web_server import run_web_simulator

__version__ = "1.0.0"
__all__ = [
    "SimulationReport",
    "TargetBindingResult",
    "BackboneAnnotationResult",
    "RCAPrimingResult",
    "ConcatemerResult",
    "QCMatrixResult",
    "find_best_arm_mapping",
    "annotate_backbone_regions",
    "simulate_rca_priming",
    "generate_rca_concatemer",
    "run_cross_reactivity_qc",
    "render_full_terminal_report",
    "generate_interactive_html",
    "run_simulation",
]
