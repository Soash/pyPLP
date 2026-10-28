"""pyPLP.multiplex.matrix_eval: Comprehensive NxN cross-reactivity interaction matrix.
"""

import csv
from pyPLP.multiplex.cross_dimer import calc_heterodimer_dg, calc_homodimer_dg


def compute_cross_reactivity_matrix(
    oligos: dict[str, str],
    max_cross_dimer_dg: float = -2.0,
    max_probe_to_noncognate_plp_dg: float = -4.0,
    max_inter_plp_dg: float = -6.0,
    temp_c: float = 37.0,
) -> dict:
    """Computes the full NxN interaction matrix across all reaction oligonucleotides

    (RCA Primer, Probes, and PLPs), evaluating thermodynamic binding free energy (dG in kcal/mol).
    Returns a dictionary containing the matrix, evaluation results, and any violations.
    """
    names = list(oligos.keys())
    matrix = {}
    violations = []

    for i, n1 in enumerate(names):
        matrix[n1] = {}
        s1 = oligos[n1]
        for j, n2 in enumerate(names):
            s2 = oligos[n2]

            if i == j:
                dg = calc_homodimer_dg(s1, temp_c=temp_c)
            else:
                dg = calc_heterodimer_dg(s1, s2, temp_c=temp_c)

            def _extract_channel(name: str) -> str:
                parts = name.split("_")
                if len(parts) >= 2 and parts[0] == "PLP" and parts[1].isdigit():
                    return parts[1]
                if len(parts) >= 3 and parts[0] == "RCA" and parts[1] == "Primer" and parts[2].isdigit():
                    return parts[2]
                if len(parts) >= 3 and parts[0] == "Detection" and parts[1] == "Probe" and parts[2].isdigit():
                    return parts[2]
                for part in parts:
                    if part.isdigit():
                        return part
                return ""

            # Determine appropriate threshold based on pair type
            is_rca_plp = ("RCA_Primer" in n1 and "PLP" in n2) or ("PLP" in n1 and "RCA_Primer" in n2)
            is_probe_plp = ("Probe" in n1 and "PLP" in n2) or ("PLP" in n1 and "Probe" in n2)

            if is_rca_plp:
                ch1 = _extract_channel(n1)
                ch2 = _extract_channel(n2)
                is_cognate = bool(ch1 and ch2 and ch1 == ch2) or (not ch1 and not ch2)
                if is_cognate:
                    threshold = -999.0
                    rule_name = "Cognate RCA Primer-to-PLP Binding"
                else:
                    threshold = max_cross_dimer_dg
                    rule_name = "Non-Cognate RCA Primer-to-PLP"
            elif is_probe_plp:
                ch1 = _extract_channel(n1)
                ch2 = _extract_channel(n2)
                is_universal = ("Universal" in n1 or "Universal" in n2)
                is_cognate = is_universal or (bool(ch1 and ch2 and ch1 == ch2))
                if is_cognate:
                    threshold = -999.0
                    rule_name = "Cognate Detection Probe"
                else:
                    threshold = max_probe_to_noncognate_plp_dg
                    rule_name = "Non-Cognate Probe-to-PLP"
            elif "PLP" in n1 and "PLP" in n2 and i != j:
                threshold = max_inter_plp_dg
                rule_name = "Inter-PLP Interaction"
            else:
                threshold = max_cross_dimer_dg
                rule_name = "Oligo Cross-Dimer"

            passes = dg >= threshold
            matrix[n1][n2] = {
                "dg": dg,
                "threshold": threshold,
                "passes": passes,
                "rule": rule_name,
            }

            if not passes and i <= j:
                violations.append({
                    "oligo1": n1,
                    "oligo2": n2,
                    "dg": dg,
                    "threshold": threshold,
                    "rule": rule_name,
                })

    return {
        "names": names,
        "matrix": matrix,
        "violations": violations,
        "is_valid": len(violations) == 0,
    }


def save_cross_reactivity_matrix_csv(matrix_data: dict, output_csv: str):
    """Exports the NxN interaction matrix as a formatted CSV table."""
    names = matrix_data["names"]
    matrix = matrix_data["matrix"]

    with open(output_csv, "w", newline="") as f_out:
        writer = csv.writer(f_out)
        # Header row: blank top-left + oligo names
        writer.writerow(["Oligo_Name"] + names)

        for n1 in names:
            row = [n1]
            for n2 in names:
                cell = matrix[n1][n2]
                val = cell["dg"]
                status = "" if cell["passes"] else " [CLASH]"
                row.append(f"{val:.2f}{status}")
            writer.writerow(row)
