/**
 * pyPLP Simulator - Presets
 */

const PRESETS = {
  mtb_dna: {
    target: "ATCGATCGATGGGATGCATGTCTTGTGGTGGAAAGCGCTTCTAGCTAGCT",
    plp: "CAAGACATGCATCCCTAAGCATAACAGCCCACCACTATGTAACCTCTAGACAACCTCCTAAGCGCTTTCCACCA",
    primer: "TAGTGGTGGGCTGTTATG",
    probe: "AACCTCTAGACAACCTCC",
    temp: 37.0,
    type: "DNA",
    title: "MTB Singleplex Assay (DNA Target)",
    desc: "High-specificity standard detection with 0-bp flush nick and balanced arm Tm."
  },
  rna_splintr: {
    target: "AUCGAUCGAUGGGAUGCAUGUCUUGUGGUGGAAAGCGCUUCUAGCUAGCU",
    plp: "CAAGACATGCATCCCTAAGCATAACAGCCCACCACTATGTAACCTCTAGACAACCTCCTAAGCGCTTTCCACCA",
    primer: "TAGTGGTGGGCTGTTATG",
    probe: "AACCTCTAGACAACCTCC",
    temp: 37.0,
    type: "RNA",
    title: "RNA Target with SplintR Ligase",
    desc: "Direct RNA detection without reverse transcription; SplintR seals DNA on RNA splint."
  },
  gap_fail: {
    target: "ATCGATCGATGGGATGCATGTCTTGATGGTGGAAAGCGCTTCTAGCTAGCT",
    plp: "CAAGACATGCATCCCTAAGCATAACAGCCCACCACTATGTAACCTCTAGACAACCTCCTAAGCGCTTTCCACCA",
    primer: "TAGTGGTGGGCTGTTATG",
    probe: "AACCTCTAGACAACCTCC",
    temp: 37.0,
    type: "DNA",
    title: "1-bp Gap (Ligation Arrest)",
    desc: "Gap between arms prevents covalent circularization; ligase cannot span gaps."
  },
  mismatch_fail: {
    target: "ATCGATCGATGGGATGCATGTCTTGCGGTGGAAAGCGCTTCTAGCTAGCT",
    plp: "CAAGACATGCATCCCTAAGCATAACAGCCCACCACTATGTAACCTCTAGACAACCTCCTAAGCGCTTTCCACCA",
    primer: "TAGTGGTGGGCTGTTATG",
    probe: "AACCTCTAGACAACCTCC",
    temp: 37.0,
    type: "DNA",
    title: "3'-End Mismatch Failure",
    desc: "Non-Watson-Crick base at the -1 position abolishes catalytic ligation efficiency."
  },
  primer_fail: {
    target: "ATCGATCGATGGGATGCATGTCTTGTGGTGGAAAGCGCTTCTAGCTAGCT",
    plp: "CAAGACATGCATCCCTAAGCATAACAGCCCACCACTATGTAACCTCTAGACAACCTCCTAAGCGCTTTCCACCA",
    primer: "CCCCCCCCCCCCCCCCCC",
    probe: "AACCTCTAGACAACCTCC",
    temp: 37.0,
    type: "DNA",
    title: "Non-Complementary RCA Primer (RCA Blocked)",
    desc: "RCA primer has no sequence complementarity to the circular padlock probe; Phi29 cannot initiate."
  },
  probe_fail: {
    target: "ATCGATCGATGGGATGCATGTCTTGTGGTGGAAAGCGCTTCTAGCTAGCT",
    plp: "CAAGACATGCATCCCTAAGCATAACAGCCCACCACTATGTAACCTCTAGACAACCTCCTAAGCGCTTTCCACCA",
    primer: "TAGTGGTGGGCTGTTATG",
    probe: "TTTTTTTTTTTTTTTTTT",
    temp: 37.0,
    type: "DNA",
    title: "Non-Complementary Detection Probe (0 a.u. Signal)",
    desc: "Detection probe has no complementary binding sites along the concatemer; probes remain free in solution."
  }
};
