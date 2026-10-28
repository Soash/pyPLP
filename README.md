# pyPLP: An Automated Computational Framework for Padlock Probe Design, Biophysical Screening, and Multiplex Isothermal Diagnostic Assays

**Author(s):** Md Abu Sadat Soash, Md. Mainul Hossain Bakul, Nahian Kabir Arpa  
**Supervisor:** [Prof. Dr. Gokul Chandra Biswas](https://www.sust.edu/departments/geb/faculty/gcbiswas-geb@sust.edu)  
*Department of Genetic Engineering and Biotechnology*<br>
*Shahjalal University of Science and Technology, Sylhet, Bangladesh*

[![Status: Manuscript in Preparation](https://img.shields.io/badge/Status-Manuscript%20in%20Preparation-orange.svg?style=flat-square)](#authorship-academic-supervision--institutional-framework)
[![Python: 3.10+](https://img.shields.io/badge/Python-3.10%2B-blue.svg?style=flat-square&logo=python)](https://www.python.org/)
[![Engines: ViennaRNA & DNA mfold](https://img.shields.io/badge/Thermodynamics-ViennaRNA%20%7C%20DNA%20mfold-purple.svg?style=flat-square)](#2-biophysical--thermodynamic-architecture)
[![Assay: RCA | LFA | CRISPR Cas12a](https://img.shields.io/badge/Assay-RCA%20%7C%20LFA%20%7C%20CRISPR%20Cas12a-brightgreen.svg?style=flat-square)](#3-crispr-cas12a-vttv-pam-engineering--shared-backbones)
[![Readout: Visual LFA & Fluorescence](https://img.shields.io/badge/Readout-Visual%20LFA%20%26%20Fluorescence-red.svg?style=flat-square)](#4-standardized-synthesis-order-sheet)

---

## Abstract and Scientific Scope

Molecular point-of-care (PoC) diagnostics require ultra-high specificity and sensitivity without dependence on complex thermal cyclers. **Padlock Probes (PLPs)** coupled with **Rolling Circle Amplification (RCA)** offer single-nucleotide discrimination and exponential isothermal signal amplification, making them an ideal foundation for instrument-free diagnostic assays such as **Lateral Flow Assays (LFA)** and **CRISPR Cas12a collateral cleavage biosensors**.

However, designing functional padlock probes presents major biophysical challenges:
1. **Target Invariance**: Pathogen genomes (e.g., Dengue virus, *Mycobacterium tuberculosis*) exhibit hypervariable regions; probes must target regions with 100% conservation across clinical isolates to avoid diagnostic escape.
2. **Unimolecular Secondary Structure**: Probes are typically 70–100 nt single-stranded oligonucleotides. Intramolecular hairpins or arm self-hybridization sequester the hybridizing arms, drastically suppressing circularization kinetics.
3. **Off-Target Cross-Reactivity**: Probes must be rigorously excluded from binding to host transcriptomes (e.g., human) or co-circulating pathogens.
4. **Multiplex Optimization**: Designing panels of probes that simultaneously avoid pairwise cross-dimers while sharing universal amplification primers and reporters is computationally demanding.

**`pyPLP`** is a high-throughput computational biology framework and automated pipeline engineered to resolve these challenges. It integrates multiple-sequence alignment (**MAFFT**), nearest-neighbor cation-dependent melting temperature thermodynamics, dual-engine secondary structure evaluation (**ViennaRNA** & **DNA mfold**), genome-wide **NCBI BLAST+** specificity screening, automated **CRISPR Cas12a VTTV PAM junction engineering**, and Monte Carlo simulated annealing optimization to generate synthesis-ready, high-performance padlock probes for infectious disease diagnostics and spatial transcriptomics.

---

## Key Computational Features

* **Dual-Engine Secondary Structure Minimization**:
  Evaluates minimum free energy ($\Delta G_{\text{MFE}}$), base-pair probabilities, and dot-bracket accessibility using both **ViennaRNA (`RNAfold`)** and **DNA mfold** (Zuker / SantaLucia parameters), guaranteeing that hybridizing arms remain accessible for target circularization.
* **Hierarchical Multi-Target Tie-Breaking Ranking**:
  Ranks candidate probes across target panels based on a composite biophysical priority matrix: ViennaRNA pass rate $\rightarrow$ ViennaRNA $\Delta G$ $\rightarrow$ mfold paired arm bases $\rightarrow$ mfold total paired bases $\rightarrow$ mfold worst $\Delta G$.
* **Automated CRISPR Cas12a (VTTV) PAM Junction Engineering**:
  Engineers canonical $[A/C/G]\text{-T-T-}[A/C/G]$ (VTTV) PAM sequences at probe junctions, assembling target-spanning crRNA guides and universal collateral reporter probes (FQ and FB reporters for fluorometry and lateral flow strips).
* **Universal Shared-Backbone Multiplex Panels**:
  Designs unified segmented backbones where diverse pathogen targets share identical RCA primer binding sites and universal detection probe sequences, radically reducing synthesis cost and multiplex complexity.
* **Genome-Wide BLAST+ Host & Cross-Reactivity Screening**:
  Automated multi-database screening against human genomic backgrounds and co-circulating clinical pathogens (e.g., CHIKV, DENV, HBV, HCV, Mtb, ZIKV).
* **In Silico Reaction Simulation**:
  Simulates padlock hybridization, enzymatic nick-ligation kinetics, and phi29 DNA polymerase concatemeric repeat yields under user-defined temperature and reaction durations.
* **Single Combined Oligo Order Export (`combined_oligos.csv`)**:
  Generates a deduplicated, synthesis-ready CSV formatted with commercial vendor modifications (`5'-Phos` for ligation, `5'-Biotin` for RCA primers, `5'-6-FAM` for fluorescent detection, and Dual-labeled FQ/FB CRISPR reporters).

---

## 1. Padlock Probe (PLP) Architecture

A padlock probe is a linear ssDNA oligonucleotide (~70–100 nt) consisting of:
* **Target-Complementary Left Arm (15–20 nt, 5'-Phos)**
* **Segmented Universal Backbone (40–60 nt)**:
  - Spacer 1 (4 nt)
  - RCA Primer Binding Region (18–20 nt)
  - Spacer 2 (4 nt)
  - Detection Probe Binding Region (18–20 nt)
  - Spacer 3 (0–4 nt)
  - Optional Standalone Restriction Site (e.g., *HaeIII*, *BstUI*)
  - Rest Segment (engineered VTTV PAM tail)
* **Target-Complementary Right Arm (15–20 nt, 3'-OH)**

```
Linear Padlock Probe (PLP) Architecture (5' -> 3'):

  5'-Phos                                                                    3'-OH
   [ Left Arm ]-[ Spacer1 ]-[ RCA Primer Site ]-[ Spacer2 ]-[ Probe Site ]-[ Rest ]-[ Right Arm ]
        |                                                                           |
        |                        Hybridization on Target ssRNA/DNA                  |
        V                                                                           V
  3'-[ Target Left ]-------------------------------------------------------------[ Target Right ]-5'
                     ^ Nick Ligation Junction (Covalently sealed by T4/SplintR Ligase)
```

Upon contiguous hybridization to the target sequence, the $5'$-phosphate and $3'$-hydroxyl termini are juxtaposed with zero nucleotide gaps, allowing a DNA/RNA ligase to circularize the probe. Subsequent Rolling Circle Amplification (RCA) catalyzed by phi29 DNA polymerase generates single-stranded tandem concatemers containing thousands of repeats of the antisense backbone.

---

## 2. Biophysical & Thermodynamic Architecture

To ensure high-efficiency hybridization and circularization, `pyPLP` filters every candidate probe through rigorous biophysical criteria:

```
+-------------------------------------------------------------------------------------------------------+
|                                    BIOPHYSICAL SCREENING CRITERIA                                     |
+------------------------------------+----------------------------------+-------------------------------+
| Parameter                          | Recommended Range / Default      | Biophysical Significance      |
+------------------------------------+----------------------------------+-------------------------------+
| Target Arm Melting Temp (Tm)       | 40.0°C – 65.0°C (ΔTm ≤ 2.0°C)   | Cation-adjusted (SantaLucia)  |
| Arm GC Content                     | 40.0% – 60.0%                    | Prevents GC/AT extreme skews  |
| Maximum Homopolymer Length         | ≤ 3–4 nt                         | Prevents synthesis slippage   |
| Unimolecular Hairpin Tm            | ≤ 30.0°C                         | Eliminates arm self-closure   |
| Maximum Arm Paired Bases           | ≤ 1–4 nt                         | Keeps hybridizing arms open   |
| Minimum Folding Energy (ΔG_MFE)    | ≥ -6.0 to -1.0 kcal/mol          | Open thermodynamic state      |
| Cross-Dimer Heterodimerization ΔG  | ≥ -5.0 kcal/mol                  | Suppresses probe-probe dimers |
| BLAST Host Alignment Identity      | ≤ 70.0% (0 allowed hits)         | Zero off-target cross-binding |
+------------------------------------+----------------------------------+-------------------------------+
```

---

## 3. CRISPR Cas12a (VTTV) PAM Engineering & Shared Backbones

In diagnostic formats combining Padlock Probes, RCA, and CRISPR Cas12a, Cas12a recognition provides a second layer of single-base specificity. Cas12a requires a **VTTV PAM** ($[A/C/G]\text{-T-T-}[A/C/G]$) positioned upstream of a 20 nt protospacer.

`pyPLP` features an automated CRISPR PAM engine:
* **Junction Engineering**: Embeds canonical VTTV PAM sequences across the circular ligation junction without altering target-specific arm sequences.
* **Universal Rest Tail (`CTTA`)**: When multiplexing targets with differing arm prefixes, `pyPLP` engineers a universal standalone VTTV PAM (`CTTA`) into the backbone Rest terminus.
* **Strict Nick-Spanning Validation**: Ensures the 20 nt crRNA spacer spans across the circular ligation nick (Right Arm + Left Arm), proving target-directed circularization while completely eliminating spurious sampling of the universal backbone.

```
Circular Ligation & Cas12a Activation Complex:

                   [ Rest: ...CTTA ] 
                          |  <- PAM (5'-CTTA-3')
                          V
  5' ... C T T A [ R1 R2 ... R15 ] | [ L1 L2 ... L5 ] L6 ... L15 [ Backbone 5' ] ... 3'
                 \_______________/   \______________/
                  15 nt Right Arm  +   5 nt Left Arm
                  ==================================
                   20 nt Target-Specific Protospacer -> crRNA Guide Activation!
```

---

## 4. Standardized Synthesis Order Sheet

`pyPLP` outputs a single, clean synthesis-ready order sheet (`combined_oligos.csv`), deduplicating universal oligos across ranks:

| Oligo Name | Modifications | Type | Functional Role |
| :--- | :--- | :--- | :--- |
| **Target_PLP_Shared_BB_1** | `5'-Phos` | DNA | Target-specific circularizing padlock probe |
| **Universal_RCA_Primer_Shared_BB_1** | `5'-Biotin` | DNA | Universal RCA amplification primer |
| **Universal_Detection_Probe_Shared_BB_1** | `5'-6-FAM` | DNA | Universal Lateral Flow / Fluorescent detection probe |
| **Target_crRNA_Shared_BB_1** | `None (RNA)` | RNA | Cas12a guide RNA (21 nt repeat + 20 nt target spacer) |
| **CRISPR_FQ_Reporter_Probe** | `5'-FAM / 3'-BHQ1` | DNA | Cas12a Fluorescent Quenched Reporter (`5'-FAM-TTATTATT-BHQ1-3'`) |
| **CRISPR_FB_Reporter_Probe** | `5'-FAM / 3'-Biotin` | DNA | Cas12a Lateral Flow Biosensor Reporter (`5'-FAM-TTATTATT-Biotin-3'`) |

---

## 5. Pipeline Architecture & Workflow

```
Input Genomic FASTA Sequences (Clinical Isolates / Target Pathogens)
                           |
                           v
               [ Step 1: MAFFT Alignment ]
             Global / Local Multi-Alignment
                           |
                           v
         [ Step 2: Conservation Extraction ]
          Pinpoint contiguous 100% conserved blocks
                           |
                           v
          [ Step 3: Candidate Arm Screening ]
        Tm balance, GC content, homopolymers, hairpins
                           |
                           v
          [ Step 4: NCBI BLAST+ Specificity ]
       Filter against Human Genome & Non-Target Organisms
                           |
                           v
      [ Step 5: Segmented Backbone Optimization ]
        Simulated Annealing Monte Carlo Generation
        ViennaRNA & DNA mfold Dual Thermodynamic Folding
                           |
                           v
      [ Step 6: Multiplex Panel & CRISPR Assembly ]
       Cross-reactivity screening, VTTV PAM engineering,
       crRNA generation, and combined oligo synthesis export
```

---

## 6. Installation & Dependencies

### Prerequisites

`pyPLP` relies on established, standard bioinformatics command-line utilities:
* **Python**: 3.10, 3.11, or 3.12
* **NCBI BLAST+**: `blastn` and `makeblastdb`
* **ViennaRNA Package**: `RNAfold`
* **UNAFold / mfold**: `mfold`, `sir_graph`
* **MAFFT**: `mafft`

### Conda Environment Setup

```bash
# Clone the repository
git clone https://github.com/Soash/pyPLP.git
cd pyPLP

# Create and activate a dedicated conda environment
conda create -n pyplp python=3.11 -y
conda activate pyplp

# Install core bioinformatics dependencies via Bioconda
conda install -c bioconda -c conda-forge mafft blast viennarna mfold biopython -y

# Install Python requirements
pip install -r requirements.txt
```

---

## 7. Quick Start & Execution

### 1. Automated End-to-End Pipeline
Execute the full pipeline from multi-FASTA input to designed multiplex probes:

```bash
python main.py \
  --input input/denv_isolates.fasta \
  --config configs/config_denv_pipeline.json \
  --output results/denv_pipeline
```

### 2. Standalone Universal Shared-Backbone Design (CRISPR Cas12a Mode)
Design universal shared backbones across a panel of pathogen target arms with integrated CRISPR Cas12a crRNAs:

```bash
python design_plp.py \
  -j denv_crispr_panel \
  -a info/DENV/denv_arms.csv \
  --crispr \
  -n 100 \
  --config info/DENV/config_denv_plp.json
```

### 3. Singleplex Padlock Probe Design
Design an optimized probe for specific target arms:

```bash
python design_plp.py \
  -j mtb_plp_design \
  --left-arm CAAGACATGCATCCC \
  --right-arm AAGCGCTTTCCACCA \
  --config configs/config_singleplex.json
```

### 4. In Silico Reaction Simulator
Simulate isothermal circularization and RCA concatemer production kinetics:

```bash
python simulate_plp.py \
  --plp-seq AGATCTGCTCTCTAATTTCCTACGTACGACACTGCCCGGCCTCTAATCATACTTTCACACTTATTGGTTGTTCATCAG \
  --target-seq CTGATGAACAACCAATTAGAGAGCAGATCT \
  --duration-minutes 60 \
  --temperature 37.0
```

---

## 8. Repository Structure

```
pyPLP/
├── pyPLP/                            # Core Python library
│   ├── alignment/                    # MAFFT alignment execution & validation
│   ├── arm_design/                   # Target arm extraction, Tm & biophysics
│   ├── blast/                        # NCBI BLAST+ background screening
│   ├── checkpoint/                   # Fault-tolerant pipeline checkpointing
│   ├── config/                       # JSON configuration parsing & schema
│   ├── conservation/                 # Shannon entropy & conserved blocks extraction
│   ├── crispr/                       # VTTV PAM junction engineering & crRNA designer
│   ├── multiplex/                    # Panel optimization & cross-dimer matrix
│   ├── pipeline/                     # End-to-end pipeline runner
│   ├── plp_generation/               # Segmented backbone optimization & structure eval
│   └── thermodynamics/               # Nearest-neighbor Tm, ViennaRNA & DNA mfold
├── configs/                          # Pre-configured assay JSON parameter profiles
├── info/                             # Target sequence references and metadata (DENV, MTB)
├── plans/                            # Permanent archival of implementation plans & walkthroughs
├── design_plp.py                     # Standalone padlock probe generator & ranking CLI
├── design_multiplex.py               # Standalone multiplex panel optimizer CLI
├── generate_universal_probe.py       # Standalone universal detection probe generator
├── simulate_plp.py                   # In silico ligation & RCA kinetics simulator
├── main.py                           # Full pipeline orchestration entry point
└── README.md                         # Project documentation
```

---

## Authorship, Academic Supervision & Institutional Framework

This computational framework and molecular diagnostic methodology is developed under the research supervision of **Prof. Dr. Gokul Chandra Biswas** at Shahjalal University of Science and Technology.

* **Author(s)**:
  - **Md Abu Sadat Soash**
  - **Md. Mainul Hossain Bakul**
  - **Nahian Kabir Arpa**
* **Research Supervisor**:
  - **[Prof. Dr. Gokul Chandra Biswas](https://www.sust.edu/departments/geb/faculty/gcbiswas-geb@sust.edu)**  
    Professor, Department of Genetic Engineering and Biotechnology (GEB)  
    Shahjalal University of Science and Technology (SUST), Sylhet–3114, Bangladesh  
    *Research Focus: Molecular Diagnostics, Infectious Disease Genomics, Biosensor Development*

---

## Citation & Manuscript Status

> **Manuscript in Preparation (2026)**  
> *pyPLP: An Automated Computational Framework for Padlock Probe Design, Biophysical Screening, and Multiplex Isothermal Diagnostic Assays.*  
> Soash, M. A. S., Bakul, M. M. H., Arpa, N. K., & Biswas, G. C.

If you utilize `pyPLP` or its design algorithms in your research, diagnostics pipeline, or academic publication, please cite this repository and framework.

---

## License

This software is developed for academic research, molecular diagnostics innovation, and educational purposes.  
Licensed under the [MIT License](LICENSE).
