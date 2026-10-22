# tests/row_digest.py

"""
A scatter's rows compared without holding them (TEST.258): since GEN.404 the toy galaxy holds a real galaxy's stellar
mass, so a shell of its scatter is hundreds of thousands to millions of rows and a test that put two of them in lists
held gigabytes. A test that compares runs, or counts them, streams the rows through these instead.
"""

import hashlib


def digest(rows, keep=None, check=None):
    """
    The number of rows and a hash of their contents, in order.

    Args:
        rows (iterable[tuple]): The rows of one scatter, drawn as they are read.
        keep (callable, optional): Only the rows it returns true for count and hash.
        check (callable, optional): Called with every row read, before `keep`, to assert something of it.

    Returns:
        tuple[int, str]: How many rows were kept and the hex digest of their `repr`.
    """
    sha = hashlib.sha256()
    count = 0
    for row in rows:
        if check is not None:
            check(row)
        if keep is None or keep(row):
            count += 1
            sha.update(repr(row).encode())
    return count, sha.hexdigest()


def count(rows):
    """How many rows there are, without keeping them or hashing them."""
    return sum(1 for _ in rows)
