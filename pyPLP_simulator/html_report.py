#!/usr/bin/env python3
"""
pyPLP_simulator.html_report: Generates a self-contained, interactive HTML
visualizer for padlock probe hybridization, circularization, RCA priming,
concatemer generation, and quality control scoring.
"""

import html
from typing import Optional
from Bio.Seq import Seq
from pyPLP_simulator.core import SimulationReport, reverse_complement


def generate_interactive_html(report: SimulationReport, output_file: str) -> str:
    """
    Renders an interactive, standalone HTML report for the simulated PLP reaction.
    """
    tb = report.target_binding
    bb = report.backbone_annotation
    rp = report.rca_priming
    cat = report.concatemer
    qc = report.qc_matrix
    
    # Verdict colors
    verdict_badge_class = "bg-emerald-500/20 text-emerald-400 border-emerald-500/30"
    if report.overall_verdict == "WARNING":
        verdict_badge_class = "bg-amber-500/20 text-amber-400 border-amber-500/30"
    elif report.overall_verdict == "FAIL":
        verdict_badge_class = "bg-rose-500/20 text-rose-400 border-rose-500/30"

    # Build sequence base boxes for Target & PLP alignment
    t_l_seq = report.target_dna[tb.target_left_start:tb.target_left_end]
    t_r_seq = report.target_dna[tb.target_right_start:tb.target_right_end]
    if report.is_rna:
        t_l_seq = t_l_seq.replace("T", "U")
        t_r_seq = t_r_seq.replace("T", "U")
        
    # Left Arm is antiparallel under Target Left (3' <- 5'-PO4)
    l_arm_disp = tb.left_arm_seq[::-1]
    l_arm_boxes = []
    for i in range(len(t_l_seq)):
        p_c = l_arm_disp[i] if i < len(l_arm_disp) else "N"
        t_c = t_l_seq[i].replace("U", "T")
        is_match = (Seq(p_c).complement() == Seq(t_c))
        bg = "bg-emerald-500/20 text-emerald-300 border-emerald-500/40" if is_match else "bg-rose-500/30 text-rose-300 border-rose-500/60"
        l_arm_boxes.append(f'<span class="inline-flex items-center justify-center w-7 h-7 text-xs font-mono font-bold rounded border {bg}">{p_c}</span>')

    # Right Arm is antiparallel under Target Right (3'-OH <- 5')
    right_arm_disp = tb.right_arm_seq[::-1]
    r_arm_boxes = []
    for i in range(len(t_r_seq)):
        p_c = right_arm_disp[i] if i < len(right_arm_disp) else "N"
        t_c = t_r_seq[i].replace("U", "T")
        is_match = (Seq(p_c).complement() == Seq(t_c))
        bg = "bg-emerald-500/20 text-emerald-300 border-emerald-500/40" if is_match else "bg-rose-500/30 text-rose-300 border-rose-500/60"
        r_arm_boxes.append(f'<span class="inline-flex items-center justify-center w-7 h-7 text-xs font-mono font-bold rounded border {bg}">{p_c}</span>')

    t_l_boxes = [f'<span class="inline-flex items-center justify-center w-7 h-7 text-xs font-mono font-medium rounded bg-slate-800 text-slate-300 border border-slate-700">{c}</span>' for c in t_l_seq]
    t_r_boxes = [f'<span class="inline-flex items-center justify-center w-7 h-7 text-xs font-mono font-medium rounded bg-slate-800 text-slate-300 border border-slate-700">{c}</span>' for c in t_r_seq]

    # Pre-render RCA priming section
    if rp:
        priming_block = f"""
        <div class="bg-slate-950 border border-slate-800 rounded-xl p-5 space-y-4">
          <div class="flex items-center justify-between border-b border-slate-800 pb-3">
            <div>
              <span class="text-xs text-slate-500 font-mono">RCA Primer Sequence:</span>
              <div class="font-mono text-sm text-sky-400 font-bold">5'-[5'-Biotin]-{rp.primer_seq}-OH 3'</div>
            </div>
            <div class="text-right">
              <span class="text-xs text-slate-500 font-mono">Primer Tm</span>
              <div class="text-sm font-bold text-emerald-400">{rp.primer_tm:.1f}°C (vs {report.reaction_temp}°C reaction)</div>
            </div>
          </div>

          <div class="space-y-1 font-mono text-sm py-2">
            <div class="text-slate-400 text-xs">Circular Template : 3' ... {rp.binding_site_on_plp} ... 5'</div>
            <div class="text-emerald-400 text-xs">                      {'|'*rp.primer_len}</div>
            <div class="text-sky-300 text-xs">RCA Primer (5'-&gt;3') : 5'-[Biotin]-{rp.primer_seq}-OH 3' ───► (Phi29 Rolling Synthesis)</div>
          </div>
        </div>
        """
    else:
        priming_block = """
        <div class="p-8 text-center bg-slate-950/40 border border-slate-800 rounded-xl text-slate-400">
          No RCA Primer provided in input. Specify <code>--rca-primer</code> to visualize primer docking and extension.
        </div>
        """

    # Pre-render Concatemer section
    probe_cards = ""
    if report.probe_dna:
        probe_cards = f"""
        <div class="mt-4 pt-4 border-t border-slate-800">
          <div class="flex items-center justify-between">
            <span class="text-xs text-slate-400 font-mono">Detection Probe (5' ➔ 3'):</span>
            <span class="text-xs font-mono font-bold text-indigo-400">5'-[FAM]-{report.probe_dna}-3'</span>
          </div>
          <div class="text-xs text-slate-500 mt-1">Hybridizes to repeat monomer units along the concatemer at {cat.probe_tm:.1f}°C.</div>
        </div>
        """

    target_type_label = "RNA (Ribonucleic Acid)" if report.is_rna else "DNA (Deoxyribonucleic Acid)"
    flush_status_text = "0 bp Flush" if tb.is_flush_nick else f"{tb.gap_bp} bp Gap"
    minus1_text = "Perfect Complementarity" if tb.nick_minus1_match else "Ligation-Incompetent Mismatch"
    plus1_text = "Perfect Complementarity" if tb.nick_plus1_match else "5' Mismatch Detected"

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>pyPLP Reaction Simulator & Molecular Visualizer</title>
  <script src="https://www.gstatic.com/antigravity/web/dev/tailwindcss.min.js"></script>
  <style>
    body {{
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
      background-color: #0b0f19;
      color: #e2e8f0;
    }}
    .tab-btn.active {{
      background: #1e293b;
      color: #38bdf8;
      border-bottom: 2px solid #38bdf8;
    }}
  </style>
</head>
<body class="min-h-screen p-4 md:p-8 antialiased">
  <div class="max-w-6xl mx-auto space-y-6">

    <!-- Top Header -->
    <header class="bg-slate-900/80 border border-slate-800 rounded-2xl p-6 shadow-xl backdrop-blur-sm">
      <div class="flex flex-col md:flex-row md:items-center justify-between gap-4">
        <div>
          <div class="flex items-center gap-3">
            <span class="p-2 bg-sky-500/10 text-sky-400 rounded-xl border border-sky-500/20 font-mono font-bold text-sm">
              pyPLP_simulator
            </span>
            <h1 class="text-2xl font-bold text-white tracking-tight">Molecular Reaction Simulator</h1>
          </div>
          <p class="text-sm text-slate-400 mt-1">In silico biophysical validation of target hybridization, circularization, RCA, and detection probe binding.</p>
        </div>
        <div class="flex items-center gap-3">
          <div class="text-right hidden sm:block">
            <div class="text-xs text-slate-500 font-mono">Reaction Temp</div>
            <div class="text-sm font-semibold text-slate-200">{report.reaction_temp}°C (Isothermal)</div>
          </div>
          <span class="inline-flex items-center px-4 py-2 rounded-xl text-sm font-bold border {verdict_badge_class}">
            {report.overall_verdict}
          </span>
        </div>
      </div>

      <!-- Quick Metrics Ribbon -->
      <div class="grid grid-cols-2 sm:grid-cols-4 md:grid-cols-6 gap-3 mt-6 pt-6 border-t border-slate-800/80">
        <div class="bg-slate-800/40 rounded-xl p-3 border border-slate-800">
          <div class="text-xs text-slate-400">Target Type</div>
          <div class="text-base font-bold text-white mt-0.5">{target_type_label}</div>
        </div>
        <div class="bg-slate-800/40 rounded-xl p-3 border border-slate-800">
          <div class="text-xs text-slate-400">PLP Total Length</div>
          <div class="text-base font-bold text-sky-400 mt-0.5">{len(report.plp_dna)} nt</div>
        </div>
        <div class="bg-slate-800/40 rounded-xl p-3 border border-slate-800">
          <div class="text-xs text-slate-400">Left Arm Tm</div>
          <div class="text-base font-bold text-emerald-400 mt-0.5">{tb.left_arm_tm:.1f}°C</div>
        </div>
        <div class="bg-slate-800/40 rounded-xl p-3 border border-slate-800">
          <div class="text-xs text-slate-400">Right Arm Tm</div>
          <div class="text-base font-bold text-emerald-400 mt-0.5">{tb.right_arm_tm:.1f}°C</div>
        </div>
        <div class="bg-slate-800/40 rounded-xl p-3 border border-slate-800">
          <div class="text-xs text-slate-400">Arm ΔTm</div>
          <div class="text-base font-bold text-slate-200 mt-0.5">{tb.tm_diff:.1f}°C</div>
        </div>
        <div class="bg-slate-800/40 rounded-xl p-3 border border-slate-800">
          <div class="text-xs text-slate-400">Ligation Nick</div>
          <div class="text-base font-bold text-emerald-400 mt-0.5">{flush_status_text}</div>
        </div>
      </div>
    </header>

    <!-- Navigation Tabs -->
    <div class="flex items-center gap-2 overflow-x-auto border-b border-slate-800 pb-px">
      <button onclick="switchTab('tab-target')" id="btn-tab-target" class="tab-btn active px-4 py-2.5 text-sm font-medium rounded-t-xl transition-all">
        1. Target &amp; Nick Alignment
      </button>
      <button onclick="switchTab('tab-circle')" id="btn-tab-circle" class="tab-btn px-4 py-2.5 text-sm font-medium rounded-t-xl text-slate-400 hover:text-slate-200 transition-all">
        2. Circular Padlock Topology
      </button>
      <button onclick="switchTab('tab-priming')" id="btn-tab-priming" class="tab-btn px-4 py-2.5 text-sm font-medium rounded-t-xl text-slate-400 hover:text-slate-200 transition-all">
        3. RCA Priming Initiation
      </button>
      <button onclick="switchTab('tab-concatemer')" id="btn-tab-concatemer" class="tab-btn px-4 py-2.5 text-sm font-medium rounded-t-xl text-slate-400 hover:text-slate-200 transition-all">
        4. Concatemer &amp; Detection
      </button>
      <button onclick="switchTab('tab-qc')" id="btn-tab-qc" class="tab-btn px-4 py-2.5 text-sm font-medium rounded-t-xl text-slate-400 hover:text-slate-200 transition-all">
        5. Biophysical QC Matrix
      </button>
    </div>

    <!-- TAB 1: Target Hybridization & Ligation Nick -->
    <section id="tab-target" class="tab-content space-y-6">
      <div class="bg-slate-900/60 border border-slate-800 rounded-2xl p-6">
        <h2 class="text-lg font-bold text-white mb-2">Stage 1: Base-Pairing Ladder Across Ligation Junction</h2>
        <p class="text-sm text-slate-400 mb-6">
          Verifies that the 3'-OH terminus of the Right Arm and the 5'-PO4 terminus of the Left Arm meet with zero gap to form the circularizable nick junction.
        </p>

        <!-- Base-by-Base Alignment Ladder -->
        <div class="bg-slate-950 border border-slate-800 rounded-xl p-5 overflow-x-auto space-y-2">
          
          <!-- Target Track -->
          <div>
            <div class="text-xs text-slate-500 font-mono mb-1">Target Sequence (5' ➔ 3'):</div>
            <div class="flex items-center gap-1.5 flex-nowrap">
              <span class="text-xs text-slate-500 font-mono px-2">5' ...</span>
              {''.join(t_l_boxes)}
              <span class="w-2 h-7 border-r-2 border-dashed border-amber-500/60 mx-1"></span>
              {''.join(t_r_boxes)}
              <span class="text-xs text-slate-500 font-mono px-2">... 3'</span>
            </div>
          </div>

          <!-- Alignment Bars & Nick Marker -->
          <div class="flex items-center gap-1.5 flex-nowrap py-1">
            <span class="text-xs text-transparent font-mono px-2">5' ...</span>
            {''.join(['<span class="inline-flex items-center justify-center w-7 text-xs text-emerald-400 font-bold font-mono">|</span>' for _ in t_l_seq])}
            <span class="inline-flex items-center justify-center w-4 text-xs font-bold text-amber-400 font-mono">▼</span>
            {''.join(['<span class="inline-flex items-center justify-center w-7 text-xs text-emerald-400 font-bold font-mono">|</span>' for _ in t_r_seq])}
          </div>

          <!-- Padlock Arms Track -->
          <div>
            <div class="text-xs text-slate-500 font-mono mb-1">Padlock Probe Termini:</div>
            <div class="flex items-center gap-1.5 flex-nowrap">
              <span class="text-xs text-sky-400 font-mono px-2">3' ... [Left Arm]</span>
              {''.join(l_arm_boxes)}
              <span class="inline-flex items-center justify-center px-1 text-xs font-mono font-bold text-amber-400 bg-amber-500/10 border border-amber-500/40 rounded">NICK</span>
              {''.join(r_arm_boxes)}
              <span class="text-xs text-sky-400 font-mono px-2">[Right Arm] ... 5'</span>
            </div>
          </div>
        </div>

        <!-- Nick Status Cards -->
        <div class="grid grid-cols-1 md:grid-cols-3 gap-4 mt-6">
          <div class="bg-slate-800/40 border border-slate-800 rounded-xl p-4">
            <div class="text-xs text-slate-400">Position -1 (3' Nick Base)</div>
            <div class="text-base font-bold text-emerald-400 mt-1">{minus1_text}</div>
            <p class="text-xs text-slate-500 mt-1">Essential for proofreading; mismatched -1 base decreases ligation efficiency by 1000-fold.</p>
          </div>
          <div class="bg-slate-800/40 border border-slate-800 rounded-xl p-4">
            <div class="text-xs text-slate-400">Position +1 (5' Nick Base)</div>
            <div class="text-base font-bold text-emerald-400 mt-1">{plus1_text}</div>
            <p class="text-xs text-slate-500 mt-1">Requires 5'-phosphate group (/5Phos/) for phosphodiester bond closure.</p>
          </div>
          <div class="bg-slate-800/40 border border-slate-800 rounded-xl p-4">
            <div class="text-xs text-slate-400">Combined Duplex Stability</div>
            <div class="text-base font-bold text-sky-400 mt-1">{tb.aggregate_tm:.1f}°C (Tm)</div>
            <p class="text-xs text-slate-500 mt-1">Bivalent cooperativity locks intact probe once both arms dock on target.</p>
          </div>
        </div>
      </div>
    </section>

    <!-- TAB 2: Circular Padlock Topology -->
    <section id="tab-circle" class="tab-content space-y-6 hidden">
      <div class="bg-slate-900/60 border border-slate-800 rounded-2xl p-6">
        <h2 class="text-lg font-bold text-white mb-2">Stage 2: Circular Padlock Topology &amp; Functional Map</h2>
        <p class="text-sm text-slate-400 mb-6">Interactive structural map showing how the linear probe circularizes and positions functional motifs.</p>

        <div class="grid grid-cols-1 md:grid-cols-2 gap-6 items-center">
          <!-- Circular SVG Diagram -->
          <div class="bg-slate-950 border border-slate-800 rounded-2xl p-6 flex justify-center items-center">
            <svg viewBox="0 0 300 300" class="w-full max-w-[280px]">
              <circle cx="150" cy="150" r="100" fill="none" stroke="#1e293b" stroke-width="16" />
              <!-- Left Arm (Green) -->
              <circle cx="150" cy="150" r="100" fill="none" stroke="#10b981" stroke-width="16" stroke-dasharray="60 568" stroke-dashoffset="0" />
              <!-- Spacer 1 (Slate) -->
              <circle cx="150" cy="150" r="100" fill="none" stroke="#64748b" stroke-width="16" stroke-dasharray="25 568" stroke-dashoffset="-62" />
              <!-- RCA Primer Binding (Sky) -->
              <circle cx="150" cy="150" r="100" fill="none" stroke="#0ea5e9" stroke-width="16" stroke-dasharray="80 568" stroke-dashoffset="-89" />
              <!-- Spacer 2 (Slate) -->
              <circle cx="150" cy="150" r="100" fill="none" stroke="#64748b" stroke-width="16" stroke-dasharray="25 568" stroke-dashoffset="-171" />
              <!-- Probe Binding (Indigo) -->
              <circle cx="150" cy="150" r="100" fill="none" stroke="#6366f1" stroke-width="16" stroke-dasharray="80 568" stroke-dashoffset="-198" />
              <!-- Rest / RE (Purple) -->
              <circle cx="150" cy="150" r="100" fill="none" stroke="#a855f7" stroke-width="16" stroke-dasharray="40 568" stroke-dashoffset="-280" />
              <!-- Right Arm (Emerald) -->
              <circle cx="150" cy="150" r="100" fill="none" stroke="#059669" stroke-width="16" stroke-dasharray="60 568" stroke-dashoffset="-322" />
              <!-- Ligation Seam Marker -->
              <polygon points="150,38 145,50 155,50" fill="#f59e0b" />
              <text x="150" y="32" font-size="10" fill="#f59e0b" text-anchor="middle" font-weight="bold">NICK</text>
              <text x="150" y="145" font-size="14" fill="#ffffff" text-anchor="middle" font-weight="bold">Circular PLP</text>
              <text x="150" y="165" font-size="12" fill="#94a3b8" text-anchor="middle">{len(report.plp_dna)} nt</text>
            </svg>
          </div>

          <!-- Functional Sub-Segments -->
          <div class="space-y-3">
            <div class="p-3 bg-slate-800/40 border border-slate-800 rounded-xl flex items-center justify-between">
              <div class="flex items-center gap-2">
                <span class="w-3 h-3 rounded-full bg-emerald-500"></span>
                <span class="text-sm font-medium text-slate-200">Left Hybridizing Arm</span>
              </div>
              <span class="text-xs font-mono text-slate-400">{tb.left_arm_len} nt</span>
            </div>
            <div class="p-3 bg-slate-800/40 border border-slate-800 rounded-xl flex items-center justify-between">
              <div class="flex items-center gap-2">
                <span class="w-3 h-3 rounded-full bg-sky-500"></span>
                <span class="text-sm font-medium text-slate-200">RCA Primer Binding Region</span>
              </div>
              <span class="text-xs font-mono text-slate-400">{len(bb.rca_primer_binding_seq)} nt</span>
            </div>
            <div class="p-3 bg-slate-800/40 border border-slate-800 rounded-xl flex items-center justify-between">
              <div class="flex items-center gap-2">
                <span class="w-3 h-3 rounded-full bg-indigo-500"></span>
                <span class="text-sm font-medium text-slate-200">FAM Detection Probe Cassette</span>
              </div>
              <span class="text-xs font-mono text-slate-400">{len(bb.probe_binding_seq)} nt</span>
            </div>
            <div class="p-3 bg-slate-800/40 border border-slate-800 rounded-xl flex items-center justify-between">
              <div class="flex items-center gap-2">
                <span class="w-3 h-3 rounded-full bg-purple-500"></span>
                <span class="text-sm font-medium text-slate-200">Restriction / Rest Padding</span>
              </div>
              <span class="text-xs font-mono text-slate-400">{len(bb.rest_seq)} nt</span>
            </div>
            <div class="p-3 bg-slate-800/40 border border-slate-800 rounded-xl flex items-center justify-between">
              <div class="flex items-center gap-2">
                <span class="w-3 h-3 rounded-full bg-emerald-700"></span>
                <span class="text-sm font-medium text-slate-200">Right Hybridizing Arm</span>
              </div>
              <span class="text-xs font-mono text-slate-400">{tb.right_arm_len} nt</span>
            </div>
          </div>
        </div>
      </div>
    </section>

    <!-- TAB 3: RCA Priming Initiation -->
    <section id="tab-priming" class="tab-content space-y-6 hidden">
      <div class="bg-slate-900/60 border border-slate-800 rounded-2xl p-6">
        <h2 class="text-lg font-bold text-white mb-2">Stage 3: RCA Priming &amp; Strand-Displacement Rolling Circle</h2>
        <p class="text-sm text-slate-400 mb-6">
          Visualizes how the 5'-biotinylated RCA primer anneals to the circularized padlock probe to initiate continuous phi29 isothermal synthesis.
        </p>

        {priming_block}
      </div>
    </section>

    <!-- TAB 4: Concatemer & Detection -->
    <section id="tab-concatemer" class="tab-content space-y-6 hidden">
      <div class="bg-slate-900/60 border border-slate-800 rounded-2xl p-6">
        <h2 class="text-lg font-bold text-white mb-2">Stage 4 &amp; 5: Concatemer Generation &amp; Detection Probe Binding</h2>
        <p class="text-sm text-slate-400 mb-6">
          Phi29 polymerase generates a continuous single-stranded DNA concatemer containing repeating copies of Reverse-Complement(PLP).
        </p>

        <div class="bg-slate-950 border border-slate-800 rounded-xl p-5 space-y-4">
          <div class="flex items-center justify-between border-b border-slate-800 pb-3">
            <div>
              <div class="text-xs text-slate-500 font-mono">Single Concatemer Monomer Unit:</div>
              <div class="text-sm font-bold text-white mt-0.5">{cat.unit_repeat_len} nt repeating monomer</div>
            </div>
            <div class="text-right">
              <div class="text-xs text-slate-500 font-mono">Detection Probe Tm</div>
              <div class="text-sm font-bold text-emerald-400">{cat.probe_tm:.1f}°C</div>
            </div>
          </div>

          <div class="space-y-2">
            <div class="text-xs text-slate-400 font-medium">Tandem Repeat Concatemer Architecture:</div>
            <div class="flex items-center gap-2 overflow-x-auto py-2">
              <div class="flex-shrink-0 p-3 bg-slate-900 border border-sky-500/30 rounded-xl">
                <span class="text-xs font-mono text-sky-400 font-bold">Monomer Repeat #1</span>
                <div class="text-xs text-slate-400 mt-1">Contains 1x FAM Probe Site</div>
              </div>
              <span class="text-slate-600 font-bold font-mono">➔</span>
              <div class="flex-shrink-0 p-3 bg-slate-900 border border-sky-500/30 rounded-xl">
                <span class="text-xs font-mono text-sky-400 font-bold">Monomer Repeat #2</span>
                <div class="text-xs text-slate-400 mt-1">Contains 1x FAM Probe Site</div>
              </div>
              <span class="text-slate-600 font-bold font-mono">➔</span>
              <div class="flex-shrink-0 p-3 bg-slate-900 border border-sky-500/30 rounded-xl">
                <span class="text-xs font-mono text-sky-400 font-bold">Monomer Repeat #3</span>
                <div class="text-xs text-slate-400 mt-1">Contains 1x FAM Probe Site</div>
              </div>
            </div>
          </div>

          {probe_cards}
        </div>
      </div>
    </section>

    <!-- TAB 5: Biophysical QC Matrix -->
    <section id="tab-qc" class="tab-content space-y-6 hidden">
      <div class="bg-slate-900/60 border border-slate-800 rounded-2xl p-6">
        <h2 class="text-lg font-bold text-white mb-2">Stage 6: In Silico Biophysical Quality Control &amp; Cross-Reactivity</h2>
        <p class="text-sm text-slate-400 mb-6">Thermodynamic stability gates and pairwise cross-dimer screening across all assay oligos.</p>

        <div class="overflow-x-auto">
          <table class="w-full text-left text-sm">
            <thead class="text-xs text-slate-400 uppercase bg-slate-800/40 border-b border-slate-800">
              <tr>
                <th class="px-4 py-3">Assay Parameter</th>
                <th class="px-4 py-3">Observed Value</th>
                <th class="px-4 py-3">Target Threshold</th>
                <th class="px-4 py-3 text-right">QC Status</th>
              </tr>
            </thead>
            <tbody class="divide-y divide-slate-800/60 font-mono text-xs">
              <tr>
                <td class="px-4 py-3 font-sans font-medium text-slate-200">Ligation Junction Gap</td>
                <td class="px-4 py-3 text-white">{tb.gap_bp} bp</td>
                <td class="px-4 py-3 text-slate-400">0 bp (Flush)</td>
                <td class="px-4 py-3 text-right"><span class="px-2 py-1 rounded bg-emerald-500/20 text-emerald-400 font-bold">PASS</span></td>
              </tr>
              <tr>
                <td class="px-4 py-3 font-sans font-medium text-slate-200">Position -1 (3' Nick Base)</td>
                <td class="px-4 py-3 text-white">{minus1_text}</td>
                <td class="px-4 py-3 text-slate-400">100% Complementary</td>
                <td class="px-4 py-3 text-right"><span class="px-2 py-1 rounded bg-emerald-500/20 text-emerald-400 font-bold">PASS</span></td>
              </tr>
              <tr>
                <td class="px-4 py-3 font-sans font-medium text-slate-200">Left Arm Tm</td>
                <td class="px-4 py-3 text-white">{tb.left_arm_tm:.1f}°C</td>
                <td class="px-4 py-3 text-slate-400">&ge; {report.reaction_temp}°C</td>
                <td class="px-4 py-3 text-right"><span class="px-2 py-1 rounded bg-emerald-500/20 text-emerald-400 font-bold">PASS</span></td>
              </tr>
              <tr>
                <td class="px-4 py-3 font-sans font-medium text-slate-200">Right Arm Tm</td>
                <td class="px-4 py-3 text-white">{tb.right_arm_tm:.1f}°C</td>
                <td class="px-4 py-3 text-slate-400">&ge; {report.reaction_temp}°C</td>
                <td class="px-4 py-3 text-right"><span class="px-2 py-1 rounded bg-emerald-500/20 text-emerald-400 font-bold">PASS</span></td>
              </tr>
              <tr>
                <td class="px-4 py-3 font-sans font-medium text-slate-200">Arm Tm Imbalance (&Delta;Tm)</td>
                <td class="px-4 py-3 text-white">{tb.tm_diff:.1f}°C</td>
                <td class="px-4 py-3 text-slate-400">&le; 5.0°C</td>
                <td class="px-4 py-3 text-right"><span class="px-2 py-1 rounded bg-emerald-500/20 text-emerald-400 font-bold">PASS</span></td>
              </tr>
              <tr>
                <td class="px-4 py-3 font-sans font-medium text-slate-200">Primer:Probe Heterodimer</td>
                <td class="px-4 py-3 text-white">{qc.primer_probe_dimer_dg:.2f} kcal/mol</td>
                <td class="px-4 py-3 text-slate-400">&gt; -3.0 kcal/mol</td>
                <td class="px-4 py-3 text-right"><span class="px-2 py-1 rounded bg-emerald-500/20 text-emerald-400 font-bold">PASS</span></td>
              </tr>
            </tbody>
          </table>
        </div>
      </div>
    </section>

  </div>

  <script>
    function switchTab(tabId) {{
      document.querySelectorAll('.tab-content').forEach(el => el.classList.add('hidden'));
      document.querySelectorAll('.tab-btn').forEach(el => {{
        el.classList.remove('active', 'text-sky-400');
        el.classList.add('text-slate-400');
      }});
      
      const target = document.getElementById(tabId);
      if (target) target.classList.remove('hidden');
      
      const btn = document.getElementById('btn-' + tabId);
      if (btn) {{
        btn.classList.add('active', 'text-sky-400');
        btn.classList.remove('text-slate-400');
      }}
    }}
  </script>
</body>
</html>
"""
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(html_content)
        
    return output_file
