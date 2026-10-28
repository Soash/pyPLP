from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pyPLP.conservation.consensus import get_column_consensus_and_conservation
from pyPLP.conservation.contiguous_blocks import find_contiguous_blocks


class ExtractionResult:
    """Result object for conserved block extraction supporting boolean evaluation,

    tuple unpacking, and property access.
    """

    def __init__(self, success: bool, actual_pct: float, block_count: int, output_file: str):
        self.success = bool(success)
        self.actual_pct = float(actual_pct)
        self.block_count = int(block_count)
        self.output_file = str(output_file)

    def __bool__(self):
        return self.success

    def __iter__(self):
        return iter((self.success, self.actual_pct, self.block_count))

    def __repr__(self):
        return (
            f"ExtractionResult(success={self.success}, actual_pct={self.actual_pct}%, "
            f"block_count={self.block_count}, output_file='{self.output_file}')"
        )


def extract_conserved_blocks(
    input_file: str,
    output_file: str,
    min_conservation_pct: float,
    min_block_len_bp: int = 20,
    ignore_gaps: bool = False,
    auto_fallback: bool = True,
    min_fallback_pct: float = 65.0,
    fallback_step: float = 1.0,
    organism_name: str = "",
) -> ExtractionResult:
    """Extracts contiguous conserved blocks from an aligned FASTA file with

    automatic percentage fallback and comprehensive metadata tracking.
    """
    alignment = AlignIO.read(input_file, "fasta")
    aln_length = alignment.get_alignment_length()
    num_seqs = len(alignment)

    org_label = organism_name or "Target"
    print(f"[Conservation] Loaded {num_seqs} sequences with alignment length {aln_length} bp for '{org_label}'.")

    # 1. Compute consensus sequence and column conservation scores
    consensus_seq = []
    scores = []
    for col_idx in range(aln_length):
        col = alignment[:, col_idx]
        base, score = get_column_consensus_and_conservation(col, ignore_gaps)
        consensus_seq.append(base)
        scores.append(score)

    # 2. Search for blocks starting at min_conservation_pct with step-by-step fallback
    curr_pct = float(min_conservation_pct)
    blocks = []
    actual_pct = curr_pct
    gap_fallback_used = False

    while curr_pct >= min_fallback_pct:
        conserved_indices = [i for i, s in enumerate(scores) if s >= curr_pct and consensus_seq[i] in "ACGT"]
        blocks = find_contiguous_blocks(conserved_indices, min_block_len_bp)
        if blocks:
            actual_pct = curr_pct
            if actual_pct < min_conservation_pct:
                print(
                    f"[Conservation] Auto-fallback triggered for '{org_label}': "
                    f"Found {len(blocks)} conserved blocks at {actual_pct:.1f}% conservation "
                    f"(requested {min_conservation_pct:.1f}%)."
                )
            break
        if not auto_fallback:
            break
        curr_pct = round(curr_pct - fallback_step, 2)

    # 3. If no blocks found and ignore_gaps was False, test gap-ignoring fallback
    if not blocks and not ignore_gaps and auto_fallback:
        print(f"[Conservation] No blocks found in strict gap mode. Testing gap-ignoring fallback for '{org_label}'...")
        scores_nogap = []
        consensus_nogap = []
        for col_idx in range(aln_length):
            b, s = get_column_consensus_and_conservation(alignment[:, col_idx], ignore_gaps=True)
            consensus_nogap.append(b)
            scores_nogap.append(s)

        curr_pct = float(min_conservation_pct)
        while curr_pct >= min_fallback_pct:
            indices = [i for i, s in enumerate(scores_nogap) if s >= curr_pct and consensus_nogap[i] in "ACGT"]
            cand_blocks = find_contiguous_blocks(indices, min_block_len_bp)
            if cand_blocks:
                blocks = cand_blocks
                scores = scores_nogap
                consensus_seq = consensus_nogap
                actual_pct = curr_pct
                gap_fallback_used = True
                print(
                    f"[Conservation] Gap-ignoring fallback successful for '{org_label}': "
                    f"Found {len(blocks)} blocks at {actual_pct:.1f}% conservation (ignore_gaps=True)."
                )
                break
            curr_pct = round(curr_pct - fallback_step, 2)

    if not blocks:
        if not auto_fallback:
            all_blocks = find_contiguous_blocks([i for i, s in enumerate(scores) if s >= min_conservation_pct and consensus_seq[i] in "ACGT"], 1)
            max_found_len = max((e - s + 1 for s, e in all_blocks), default=0)
            print(
                f"\n[Conservation Filters Applied] Strict {min_conservation_pct:.1f}% conservation, min length >= {min_block_len_bp} bp (Auto-Fallback: DISABLED, ignore_gaps={ignore_gaps}).\n"
                f"[Conservation] No contiguous blocks found for '{org_label}' meeting both criteria.\n"
                f"  -> Diagnostics: Longest contiguous block found at {min_conservation_pct:.1f}% was {max_found_len} bp (required >= {min_block_len_bp} bp).\n"
                f"  -> Solution: Set 'auto_fallback': true in your config or pass '-p <pct>' with a lower conservation threshold (e.g. -p 75.0)."
            )
        else:
            print(
                f"\n[Conservation Filters Applied] Tested conservation from {min_conservation_pct:.1f}% down to {min_fallback_pct:.1f}% (step: {fallback_step}%), min block length >= {min_block_len_bp} bp (ignore_gaps={ignore_gaps}).\n"
                f"[Conservation] No blocks found for '{org_label}' down to {min_fallback_pct:.1f}% conservation and min length >= {min_block_len_bp} bp."
            )
        return ExtractionResult(False, 0.0, 0, output_file)

    output_records = []
    for i, (start_idx, end_idx) in enumerate(blocks, start=1):
        block_len = end_idx - start_idx + 1
        subseq = "".join(consensus_seq[start_idx : end_idx + 1])
        block_scores = scores[start_idx : end_idx + 1]
        avg_score = sum(block_scores) / block_len
        min_score = min(block_scores)

        pos_str = f"aln_pos={start_idx + 1}-{end_idx + 1}"
        len_str = f"length={block_len}"
        cons_str = f"conservation_pct={actual_pct:.1f} avg_conservation_pct={avg_score:.2f} min_conservation_pct={min_score:.2f}"
        org_str = f"organism={organism_name}" if organism_name else ""
        gap_str = "gap_mode=ignore_gaps" if gap_fallback_used or ignore_gaps else "gap_mode=strict"

        rec = SeqRecord(
            Seq(subseq),
            id=f"conserved_region_{i}",
            description=f"{pos_str} {len_str} {org_str} {cons_str} {gap_str}".strip(),
        )
        output_records.append(rec)

    total_blocks = len(output_records)
    lengths = [len(r.seq) for r in output_records]
    total_bases = sum(lengths)
    coverage_pct = (total_bases / aln_length) * 100.0

    print(f"\n--- Alignment Conservation Metrics [{org_label}] ---")
    print(f"Conservation Threshold  : {actual_pct:.1f}% {'(Auto-Fallback)' if actual_pct < min_conservation_pct else ''}")
    print(f"Total Conserved Blocks  : {total_blocks}")
    print(f"Total Conserved Bases   : {total_bases} bp")
    print(f"Alignment Coverage      : {coverage_pct:.2f}% of alignment")
    print(f"Maximum Block Length    : {max(lengths)} bp")
    print(f"Minimum Block Length    : {min(lengths)} bp")
    print(f"Average Block Length    : {total_bases / total_blocks:.2f} bp")
    print("---------------------------------------------------\n")

    SeqIO.write(output_records, output_file, "fasta")
    print(f"[Conservation] Saved position-ordered blocks to: {output_file}")

    sorted_output_file = output_file.replace(".fasta", "_sorted_by_size.fasta")
    sorted_records = sorted(output_records, key=lambda x: len(x.seq), reverse=True)
    SeqIO.write(sorted_records, sorted_output_file, "fasta")
    print(f"[Conservation] Saved size-sorted blocks to: {sorted_output_file}")

    return ExtractionResult(True, actual_pct, total_blocks, output_file)


# Alias for backwards compatibility with extract_conserved.py
extract_blocks = extract_conserved_blocks
