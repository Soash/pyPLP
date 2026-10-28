#!/usr/bin/env python3
"""Main entry point for pyPLP (Padlock Probe Design Pipeline).

Usage Examples:
---------------
1. Run the entire pipeline via CLI:
   python main.py -j test_job_1 -f input/sequences.fasta -c config.json

2. Run with forced re-execution (ignoring previous checkpoints):
   python main.py -j test_job_1 -f input/sequences.fasta -c config.json --force

3. Programmatic usage in Python:
   import pyPLP

   # Run the full pipeline
   pyPLP.run_pipeline(
       job_name="test_job_1",
       fasta_file="input/sequences.fasta",
       config_file="config.json"
   )

   # Or access individual modules directly:
   from pyPLP.thermodynamics import compute_metrics, calculate_tm_oligoanalyzer
   from pyPLP.conservation import extract_conserved_blocks
   from pyPLP.arm_design import design_plp_arms
   from pyPLP.blast import discover_blast_dbs
"""

import sys
from pyPLP.pipeline.runner import main, run_pipeline

if __name__ == "__main__":
    main()
