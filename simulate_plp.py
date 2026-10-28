#!/usr/bin/env python3
"""
simulate_plp.py: Top-level entrypoint for pyPLP_simulator.
Simulates and visualizes padlock probe target hybridization, ligation, circularization,
RCA priming, concatemer generation, and detection probe binding in silico.
"""

import os
import sys

# Ensure PLP root directory is in sys.path
root_dir = os.path.dirname(os.path.abspath(__file__))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from pyPLP_simulator.cli import main

if __name__ == "__main__":
    main()
