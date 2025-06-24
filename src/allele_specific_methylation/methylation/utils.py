from typing import Iterable


def hamming(s0: Iterable, s1: Iterable) -> int:
    """Find the hamming distance between two sequences

    :param s0: first sequence
    :param s1: second sequence
    :return: the hamming distance between the two sequences
    """
    return sum(c0 != c1 for c0, c1 in zip(s0, s1, strict=False))
