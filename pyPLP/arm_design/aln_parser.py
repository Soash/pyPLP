import re


def parse_aln_pos_offset(header_description: str):
    """Extracts (aln_start, aln_end) from headers like:

    >conserved_region_92 aln_pos=11572-11617 length=46
    """
    match = re.search(r"aln_pos=(\d+)-(\d+)", header_description)
    if match:
        return int(match.group(1)), int(match.group(2))
    return 1, None


def parse_conservation_metadata(header_description: str) -> dict:
    """Extracts conservation percentage and organism metadata from FASTA headers."""
    match_cons = re.search(r"conservation_pct=([\d\.]+)", header_description)
    match_avg = re.search(r"avg_conservation_pct=([\d\.]+)", header_description)
    match_org = re.search(r"organism=([^\s]+)", header_description)
    return {
        "conservation_pct": float(match_cons.group(1)) if match_cons else None,
        "avg_conservation_pct": float(match_avg.group(1)) if match_avg else None,
        "organism": match_org.group(1) if match_org else "",
    }

