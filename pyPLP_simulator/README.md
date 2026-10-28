# pyPLP_simulator: In Silico Padlock Probe (PLP) & RCA Reaction Simulator

`pyPLP_simulator` is a dedicated molecular assay simulator and biophysical validator designed to test, verify, and visualize how a target sequence, padlock probe (PLP), RCA amplification primer, and detection probe will interact *in silico*.

---

## Key Features

1. **Target Hybridization & Ligation Nick Validation**:
   - Maps Left and Right arms onto Target sequence (supports both DNA and RNA).
   - Verifies 0-bp flush nick juxtaposition ($3'\text{-OH}$ of Right Arm directly flush against $5'\text{-PO}_4$ of Left Arm).
   - Evaluates terminal fidelity: proofreading $-1$ base check and $+1$ base check.
   - Computes individual arm melting temperatures ($T_m$), free energies ($\Delta G^\circ$), and aggregate duplex stability.
2. **Circular Padlock Probe Topology**:
   - Circularizes probe in silico and annotates functional sub-segments (Spacers, RCA binding site, FAM probe cassette, Restriction endonuclease site, and Rest padding).
   - Generates interactive circular SVG maps.
3. **RCA Priming & Strand Displacement**:
   - Verifies RCA primer annealing orientation ($3'$ facing the replication direction).
   - Simulates rolling circle isothermal synthesis around circular PLP template.
4. **Concatemer Generation & Detection Probe Hybridization**:
   - Synthesizes theoretical tandem single-stranded repeat concatemers: $[\text{Reverse Complement}(\text{PLP})]_N$.
   - Quantifies repeating detection probe (e.g. 5'-FAM) hybridization frequency and stability.
5. **Cross-Reactivity QC Matrix**:
   - Screens pairwise heterodimers (Primer:Probe, PLP:Primer, PLP:Probe).
   - Evaluates unimolecular hairpins and homodimers.
   - Issues an automated **PASS / WARNING / FAIL** quality verdict.
6. **Dual Visualizations**:
   - ASCII base-pairing ladders directly in the terminal.
   - Self-contained, responsive, tabbed **interactive HTML visualizer**.

---

## Quickstart & CLI Usage

### Direct Command Line
```bash
python simulate_plp.py \
  --target "ATCGATCGATGGGATGCATGTCTTGTGGTGGAAAGCGCTTCTAGCTAGCT" \
  --plp "CAAGACATGCATCCCTAAACATAACAGCCCACCACTACACCAACCTCTAGACAACCTCCCAGCCACCATACACTTAAGCGCTTTCCACCA" \
  --rca-primer "TAGTGGTGGGCTGTTATG" \
  --probe "AACCTCTAGACAACCTCC" \
  -T 37.0 \
  -j my_plp_simulation
```

### RNA Target
```bash
python simulate_plp.py \
  --target "AUCGAUCGAUGGGAUGCAUGUCUUGUGGUGGAAAGCGCUUCUAGCUAGCU" \
  --plp "CAAGACATGCATCCCTAAACATAACAGCCCACCACTACACCAACCTCTAGACAACCTCCCAGCCACCATACACTTAAGCGCTTTCCACCA" \
  --rca-primer "TAGTGGTGGGCTGTTATG" \
  --probe "AACCTCTAGACAACCTCC" \
  --rna \
  -T 37.0
```

### Interactive Mode
Run without arguments (or with `-i / --interactive`) to enter sequences step-by-step:
```bash
python simulate_plp.py --interactive
```

### Animated Web App Simulator
Launch the interactive 60fps molecular reaction studio with live animations of target hybridization, flush-nick ligation, continuous phi29 rolling circle amplification, and fluorometer curves:
```bash
python simulate_plp.py --web
# Or specify a custom port:
python simulate_plp.py --web --port 8080
```
Or open `pyPLP_simulator/web/web_simulator.html` directly in any web browser offline.

### Python API
```python
from pyPLP_simulator import run_simulation

report = run_simulation(
    target="ATCGATCGATGGGATGCATGTCTTGTGGTGGAAAGCGCTTCTAGCTAGCT",
    plp="CAAGACATGCATCCCTAAACATAACAGCCCACCACTACACCAACCTCTAGACAACCTCCCAGCCACCATACACTTAAGCGCTTTCCACCA",
    rca_primer="TAGTGGTGGGCTGTTATG",
    probe="AACCTCTAGACAACCTCC",
    reaction_temp=37.0,
    job_name="my_simulation",
)

print(report.overall_verdict)
print(f"Left Arm Tm: {report.target_binding.left_arm_tm}°C")
print(f"Right Arm Tm: {report.target_binding.right_arm_tm}°C")
```

---

## CLI Options

| Option | Shorthand | Description | Default |
| :--- | :---: | :--- | :--- |
| `--target` | `-t` | Target sequence ($5' \rightarrow 3'$, DNA/RNA) or path to FASTA file. | *Required* |
| `--plp` | `-p` | Padlock Probe sequence ($5' \rightarrow 3'$) or path to file. | *Required* |
| `--rca-primer` | `-r` | RCA amplification primer sequence ($5' \rightarrow 3'$). | `None` |
| `--probe` | `-d` | Detection probe sequence ($5' \rightarrow 3'$, e.g. 5'-FAM labeled). | `None` |
| `--temp` | `-T` | Isothermal reaction temperature in °C. | `37.0` |
| `--rna` | | Explicitly declare target is RNA (converts `U` to `T` internally while displaying `U`). | Auto-detected |
| `--job-name` | `-j` | Job name / directory under `results/`. | `.` |
| `--output-html` | `-o` | Custom HTML report path. | `simulation_report.html` |
| `--json` | | Path to export machine-readable JSON summary report. | `None` |
| `--quiet` | `-q` | Suppress terminal printout; generate report files only. | `False` |
| `--interactive` | `-i` | Launch interactive sequence entry prompt. | `False` |

---

## Canvas Text & Annotation Configuration

You can customize, reword, or localize all textual content displayed across the canvas (such as strand labels, nick badges, enzyme names, failure alerts, and stage descriptions) using two convenient methods:

### Method 1: Edit `canvas_config.js` Directly
All canvas text definitions are organized into clear categories in `pyPLP_simulator/web/canvas_config.js` (also provided as `canvas_config.json`):

```javascript
const CANVAS_TEXT_CONFIG = {
  strands: {
    targetStrandLabel: "Target ({targetType}) 5' → 3'",
    leftArmLinearLabel: "5'-Left Arm",
    rightArmLinearLabel: "Right Arm-3'",
    backboneLinearLabel: "Backbone",
    leftArmHybridLabel: "Left Arm (5'-PO4)",
    rightArmHybridLabel: "Right Arm (3'-OH)",
    backboneLoopLabel: "Backbone Loop",
    nucleotide5PrimeTag: "5'",
    nucleotide3PrimeTag: "3'"
  },
  ligation: {
    nickFlushBadge: "0-bp Nick (3'-OH | 5'-PO4)",
    nickMismatchBadge: "3'-End Mismatch (Arrested)",
    nickGapBadge: "{gap}-bp Gap (Arrested)",
    nickNoAlignBadge: "No Alignment",
    enzymeDnaLigase: "Ligase",
    enzymeRnaLigase: "SplintR",
    sealedSuccessBanner: "COVALENT PHOSPHODIESTER BOND SEALED",
    // ...
  },
  rcaPrimer: {
    circularPlpTitle: "Circular Padlock Probe (Closed Single Strand)",
    primerAnnealedBanner: "RCA PRIMER ANNEALED (3'-OH EXTENSION READY)",
    primerStrandLabel: "RCA Primer (5' → 3')",
    primerBlockedBadge: "✕ NON-COMPLEMENTARY PRIMER",
    // ...
  },
  rollingCircle: {
    phi29EnzymeName: "phi29",
    concatemerStrandLabel: "Single-Stranded Concatemer DNA (Antisense)",
    repeatBadgePrefix: "Repeat {repeatNum}",
    // ...
  },
  detection: {
    fluorophoreActiveBadge: "FAM",
    fluorophoreQuenchedBadge: "FAM (Off)",
    probeStrandLabel: "Detection Probe (5' → 3')",
    fluorometerIntensityUnit: "a.u.",
    // ...
  },
  stages: {
    1: { title: "Stage 1: Solution State", desc: "..." },
    2: { title: "Stage 2: Target Recognition & Hybridization", desc: "..." },
    // ...
  }
};
```

#### Dynamic Placeholders
The configuration supports dynamic placeholders that are populated during live simulation:
- `{targetType}`: Molecule type (`DNA` or `RNA`)
- `{gap}`: Ligation gap distance in base pairs (e.g., `1`, `2`)
- `{matchPct}`: Sequence alignment match percentage (e.g., `94`)
- `{repeatNum}`: Concatemer repeat index (e.g., `1`, `2`, `3`)
- `{reason}`: Diagnostic arrest reason string

### Method 2: In-Browser GUI Manager
1. In `web_simulator.html`, open the **Settings Drawer** (gear icon in the top header).
2. Click **Adjust Canvas Text Content** under *Canvas Annotations & Text Content*.
3. Edit any label in the tabbed interface (Strands, Ligation, RCA Primer, phi29, Detection, Stages).
4. Click **Apply Changes Live** to update the canvas instantly without reloading.
5. Click **Export canvas_config.js** to download the customized file to your computer.

