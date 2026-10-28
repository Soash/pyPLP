"""pyPLP: A modular Python package for automated padlock probe design and screening.

"""

__version__ = "0.1.0"

import os
import sys

# Ensure active python binary dir and miniconda bin are in PATH for mfold / sir_graph / RNAfold
conda_bin = os.path.expanduser("~/miniconda3/bin")
prefix_bin = os.path.join(sys.prefix, "bin")
current_path = os.environ.get("PATH", "")
paths_to_add = [p for p in [prefix_bin, conda_bin] if os.path.isdir(p) and p not in current_path]
if paths_to_add:
    os.environ["PATH"] = os.pathsep.join(paths_to_add) + os.pathsep + current_path

from pyPLP.checkpoint.manager import CheckpointManager
from pyPLP.alignment.mafft import run_mafft
from pyPLP.conservation.extract_blocks import extract_conserved_blocks, extract_blocks
from pyPLP.conservation.consensus import get_column_consensus_and_conservation
from pyPLP.conservation.contiguous_blocks import find_contiguous_blocks
from pyPLP.arm_design.designer import design_plp_arms
from pyPLP.config.loader import load_config
from pyPLP.pipeline.runner import run_pipeline

from pyPLP.multiplex.panel_optimizer import design_multiplex_panel

from pyPLP import alignment
from pyPLP import arm_design
from pyPLP import blast
from pyPLP import checkpoint
from pyPLP import config
from pyPLP import conservation
from pyPLP import multiplex
from pyPLP import pipeline
from pyPLP import plp_generation
from pyPLP import thermodynamics

__all__ = [
    "__version__",
    "CheckpointManager",
    "run_mafft",
    "extract_conserved_blocks",
    "extract_blocks",
    "get_column_consensus_and_conservation",
    "find_contiguous_blocks",
    "design_plp_arms",
    "design_multiplex_panel",
    "load_config",
    "run_pipeline",
    "alignment",
    "arm_design",
    "blast",
    "checkpoint",
    "config",
    "conservation",
    "multiplex",
    "pipeline",
    "plp_generation",
    "thermodynamics",
]
