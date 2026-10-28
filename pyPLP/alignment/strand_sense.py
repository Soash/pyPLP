"""pyPLP.alignment.strand_sense: Preprocessing of genomic strand sense (vRNA vs cRNA)
for negative-sense ssRNA viruses and strand-specific padlock probe targeting.
"""

import os
from Bio import SeqIO
from Bio.Seq import Seq


def resolve_target_strand_fasta(
    input_fasta: str,
    target_strand: str = "positive",
    output_fasta: str | None = None,
) -> str:
    """Preprocesses input FASTA sequences based on target strand sense.

    Args:
        input_fasta: Path to input raw FASTA file.
        target_strand: Target sense convention:
            - 'positive': Standard sense as provided in FASTA (positive-sense ssRNA, mRNA, or dsDNA).
            - 'negative': Negative-sense viral genomic RNA (vRNA). Automatically generates
                          reverse-complement sequences (antigenomic cRNA / mRNA) so that
                          padlock probes capture the abundant viral transcriptomic RNA.
            - 'both': Generates both vRNA and cRNA sequences for dual-strand evaluation.
        output_fasta: Optional destination path. If not provided, an appropriately named
                      file in the same directory or temporary location will be used.

    Returns:
        Path to the resolved FASTA file to be passed downstream to alignment and screening.
    """
    mode = str(target_strand).strip().lower()
    if mode in ("positive", "pos", "sense", "+", "vrna_as_is"):
        return input_fasta

    if not os.path.isfile(input_fasta):
        raise FileNotFoundError(f"Input FASTA file not found: '{input_fasta}'")

    if not output_fasta:
        base_dir = os.path.dirname(os.path.abspath(input_fasta))
        stem = os.path.splitext(os.path.basename(input_fasta))[0]
        if mode in ("negative", "neg", "crna", "-"):
            output_fasta = os.path.join(base_dir, f"{stem}_cRNA.fasta")
        else:
            output_fasta = os.path.join(base_dir, f"{stem}_dual_strand.fasta")

    records = list(SeqIO.parse(input_fasta, "fasta"))
    if not records:
        raise ValueError(f"No FASTA records found in '{input_fasta}'.")

    processed_records = []
    if mode in ("negative", "neg", "crna", "-"):
        for rec in records:
            rc_seq = rec.seq.reverse_complement()
            rec_id = f"{rec.id}_cRNA"
            rec_desc = f"{rec.description} [cRNA / antigenomic mRNA complement]"
            new_rec = rec.__class__(rc_seq, id=rec_id, description=rec_desc)
            processed_records.append(new_rec)
    elif mode in ("both", "dual", "all"):
        for rec in records:
            # Add original vRNA
            vrna_id = f"{rec.id}_vRNA"
            vrna_rec = rec.__class__(rec.seq, id=vrna_id, description=f"{rec.description} [vRNA]")
            processed_records.append(vrna_rec)

            # Add reverse-complement cRNA
            rc_seq = rec.seq.reverse_complement()
            crna_id = f"{rec.id}_cRNA"
            crna_rec = rec.__class__(rc_seq, id=crna_id, description=f"{rec.description} [cRNA]")
            processed_records.append(crna_rec)
    else:
        # Fallback to positive if unknown
        return input_fasta

    os.makedirs(os.path.dirname(os.path.abspath(output_fasta)), exist_ok=True)
    with open(output_fasta, "w", encoding="utf-8") as f_out:
        SeqIO.write(processed_records, f_out, "fasta")

    return output_fasta
