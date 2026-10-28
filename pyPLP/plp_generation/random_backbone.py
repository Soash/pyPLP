import random


def generate_random_backbone(length: int) -> str:
    """Generates a random nucleotide sequence of specified length."""
    return "".join(random.choices(["A", "C", "G", "T"], k=length))
