"""pyPLP.multiplex: Orthogonal multiplex padlock probe and oligo panel design.
"""

from pyPLP.multiplex.cross_dimer import (
    calc_heterodimer_dg,
    calc_homodimer_dg,
    eval_oligo_orthogonality,
)
from pyPLP.multiplex.probe_generator import (
    generate_orthogonal_probe_set,
    generate_orthogonal_rca_primers,
    generate_universal_detection_probe,
)
from pyPLP.multiplex.matrix_eval import compute_cross_reactivity_matrix
from pyPLP.multiplex.panel_optimizer import design_multiplex_panel
from pyPLP.multiplex.fixed_panel import parse_fixed_oligos, evaluate_fixed_channel

__all__ = [
    "calc_heterodimer_dg",
    "calc_homodimer_dg",
    "eval_oligo_orthogonality",
    "generate_orthogonal_probe_set",
    "generate_orthogonal_rca_primers",
    "generate_universal_detection_probe",
    "compute_cross_reactivity_matrix",
    "design_multiplex_panel",
    "parse_fixed_oligos",
    "evaluate_fixed_channel",
]
