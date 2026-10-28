import subprocess
from pyPLP.blast.junction_eval import evaluate_hit_ligation_competence
from pyPLP.blast.word_size import get_optimal_word_size


def run_blastn_short_parallel(
    unique_queries: dict,
    db_prefix: str,
    arm_len: int,
    min_identity: float = 90.0,
    threads: int = 32,
    check_ligation_junction: bool = False,
) -> dict | tuple[dict, dict]:
    """Batches unique oligonucleotide queries across all CPU threads using blastn-short.

    unique_queries: dict of {unique_seq_str: query_id}
    Returns:
        dict of {unique_seq_str: hit_count} if check_ligation_junction is False,
        or tuple of ({unique_seq_str: hit_count}, {unique_seq_str: lig_competent_hit_count})
        if check_ligation_junction is True.
    """
    if not unique_queries:
        if check_ligation_junction:
            return {}, {}
        return {}

    fasta_payload = "\n".join(
        [f">{qid}\n{seq}" for seq, qid in unique_queries.items()]
    )
    word_size = get_optimal_word_size(arm_len)
    outfmt = "6 qseqid length pident qstart qend qseq sseq" if check_ligation_junction else "6 qseqid length pident"

    cmd = [
        "blastn",
        "-task",
        "blastn-short",
        "-db",
        db_prefix,
        "-num_threads",
        str(threads),
        "-word_size",
        str(word_size),
        "-evalue",
        "1000",
        "-dust",
        "no",
        "-max_target_seqs",
        "500",  # Protect against memory explosion in high-repetition zones
        "-outfmt",
        outfmt,
    ]

    result = subprocess.run(
        cmd,
        input=fasta_payload,
        text=True,
        capture_output=True,
        check=False,
    )

    qid_to_seq = {qid: seq for seq, qid in unique_queries.items()}
    hit_counts = {seq: 0 for seq in unique_queries}
    lig_counts = {seq: 0 for seq in unique_queries}

    if result.returncode != 0:
        print(
            f"Warning running BLAST against {db_prefix}: {result.stderr.strip()[:150]}"
        )
        if check_ligation_junction:
            return hit_counts, lig_counts
        return hit_counts

    # Near full-length off-target threshold: allow max 1 nt truncation/mismatch
    min_length_match = arm_len - 1
    single_arm_len = arm_len // 2

    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        qid = parts[0]
        match_len = int(parts[1])
        pident = float(parts[2])

        if match_len >= min_length_match and pident >= min_identity:
            seq = qid_to_seq.get(qid)
            if seq:
                hit_counts[seq] += 1
                if check_ligation_junction and len(parts) >= 7:
                    qstart = int(parts[3])
                    qend = int(parts[4])
                    qseq = parts[5]
                    sseq = parts[6]
                    is_competent, _ = evaluate_hit_ligation_competence(
                        qstart, qend, qseq, sseq, arm_len=single_arm_len
                    )
                    if is_competent:
                        lig_counts[seq] += 1

    if check_ligation_junction:
        return hit_counts, lig_counts
    return hit_counts
