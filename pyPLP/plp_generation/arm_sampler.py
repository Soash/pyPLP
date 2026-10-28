import csv
import os
import random


def sample_candidate_arms(
    candidate_csv: str,
    sample_count: int = 5,
    random_seed: int | None = None,
) -> list[dict]:
    """Reads candidate_plp_arms.csv and randomly samples N candidate arm pairs.

    If candidate count is less than or equal to sample_count, all candidates are returned.
    """
    if not os.path.isfile(candidate_csv):
        raise FileNotFoundError(f"Candidate arms CSV not found at '{candidate_csv}'")

    candidates = []
    with open(candidate_csv, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            candidates.append(row)

    if not candidates:
        return []

    if random_seed is not None:
        rng = random.Random(random_seed)
    else:
        rng = random.Random()

    if len(candidates) <= sample_count:
        return candidates

    return rng.sample(candidates, sample_count)
