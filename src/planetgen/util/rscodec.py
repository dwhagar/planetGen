# planetgen/util/rscodec.py

"""
Reed-Solomon Codec
==================

A Reed-Solomon erasure code over GF(256) in the standard library (DB.43, the
codec of the sector parity file, docs/design/db-check-and-parity-repair.md
section 5.9). `RS(k, m)` makes `m` parity blocks from `k` data blocks, and any
`k` of the `k + m` blocks give the data back, so up to `m` lost blocks are
rebuilt exactly. `k + m` may not pass 256.

The generator is a Cauchy matrix scaled so row 0 is all ones: the first parity
block is the plain XOR of the data blocks, so `m = 1` is ordinary XOR parity.
Every block operation works on whole blocks as Python integers (a byte-table
translate to multiply by a constant, then one big XOR), so a block costs a
handful of C-speed calls, not a loop over its bytes.

Blocks of different lengths work through zero padding: `encode` pads every
block to the longest, the parity blocks have that length, and `decode` returns
blocks of that length for the caller to trim to the lengths it stored. `update`
makes the new parity when one member changes, from the old and new bytes alone.

The multiplication tables and the Cauchy matrices are built once per process
and shared.
"""

import functools

_POLYNOMIAL = 0x11D
"""int: The field's reducing polynomial, x^8 + x^4 + x^3 + x^2 + 1."""

MAX_BLOCKS = 256
"""int: The most blocks, data and parity together, a code can have."""

_EXP = [0] * 512
_LOG = [0] * 256


def _build_logs():
    value = 1
    for power in range(255):
        _EXP[power] = value
        _LOG[value] = power
        value <<= 1
        if value & 0x100:
            value ^= _POLYNOMIAL
    for power in range(255, 512):
        _EXP[power] = _EXP[power - 255]


_build_logs()


def _multiply(a, b):
    return 0 if a == 0 or b == 0 else _EXP[_LOG[a] + _LOG[b]]


def _inverse(a):
    return _EXP[255 - _LOG[a]]


@functools.lru_cache(maxsize=None)
def _table(constant):
    """The 256-byte table that multiplies every byte by `constant`, for `bytes.translate`."""
    return bytes(_multiply(constant, value) for value in range(256))


@functools.lru_cache(maxsize=None)
def _cauchy(k, m):
    """The `m` x `k` Cauchy matrix with row 0 scaled to all ones, as a tuple of tuples."""
    matrix = [[_inverse((k + j) ^ i) for i in range(k)] for j in range(m)]
    for i in range(k):
        scale = _inverse(matrix[0][i])
        for j in range(m):
            matrix[j][i] = _multiply(matrix[j][i], scale)
    return tuple(tuple(row) for row in matrix)


def _invert(matrix):
    """The inverse of a square matrix over GF(256) by Gauss-Jordan elimination; `ValueError` if it is singular."""
    n = len(matrix)
    work = [list(row) + [int(i == j) for j in range(n)] for i, row in enumerate(matrix)]
    for col in range(n):
        pivot = next((r for r in range(col, n) if work[r][col]), None)
        if pivot is None:
            raise ValueError("the blocks do not determine the data (a singular matrix)")
        work[col], work[pivot] = work[pivot], work[col]
        scale = _inverse(work[col][col])
        work[col] = [_multiply(scale, value) for value in work[col]]
        for row in range(n):
            if row != col and work[row][col]:
                factor = work[row][col]
                work[row] = [work[row][i] ^ _multiply(factor, work[col][i]) for i in range(2 * n)]
    return [row[n:] for row in work]


def _scaled(block, constant):
    """`block` with every byte multiplied by `constant`, as an integer."""
    if constant == 1:
        return int.from_bytes(block, "little")
    return int.from_bytes(block.translate(_table(constant)), "little")


def _combine(coefficients, blocks, length):
    """The sum over GF(256) of `coefficient * block`, as `length` bytes (shorter blocks count as zero-padded)."""
    total = 0
    for coefficient, block in zip(coefficients, blocks):
        if coefficient:
            total ^= _scaled(block, coefficient)
    return total.to_bytes(length, "little")


class RS:
    """
    An erasure code of `k` data blocks and `m` parity blocks.

    Parity block `j` has index `k + j` in `decode`; data block `i` has index `i`.

    Args:
        k (int): Data blocks, at least 1.
        m (int): Parity blocks, at least 1; `k + m` is at most `MAX_BLOCKS`.
    """

    def __init__(self, k, m):
        if k < 1 or m < 1 or k + m > MAX_BLOCKS:
            raise ValueError(f"RS needs k >= 1, m >= 1 and k + m <= {MAX_BLOCKS}, not k={k}, m={m}")
        self.k = k
        self.m = m
        self._matrix = _cauchy(k, m)

    def __repr__(self):
        return f"RS(k={self.k}, m={self.m})"

    def encode(self, blocks):
        """
        The parity blocks of a group.

        Args:
            blocks (list of bytes): The `k` data blocks; shorter ones count as zero-padded to the longest.

        Returns:
            list of bytes: `m` parity blocks, each as long as the longest data block. Block 0 is the XOR of the data.
        """
        if len(blocks) != self.k:
            raise ValueError(f"encode takes {self.k} blocks, not {len(blocks)}")
        length = max(len(block) for block in blocks)
        return [_combine(self._matrix[j], blocks, length) for j in range(self.m)]

    def _row(self, index):
        """The row of the code's generator matrix that block `index` is: an identity row for data, a Cauchy row for parity."""
        if index < self.k:
            return [int(column == index) for column in range(self.k)]
        return list(self._matrix[index - self.k])

    def decode(self, have):
        """
        The data blocks, from any `k` blocks of the group.

        Args:
            have (dict): Index to bytes (data blocks `0..k-1`, parity blocks `k..k+m-1`) for at least `k` blocks that
                survived. A block shorter than the others counts as zero-padded.

        Returns:
            list of bytes: The `k` data blocks, each as long as the longest block given (trim each to the length you
                stored for it).
        """
        if len(have) < self.k:
            raise ValueError(f"decode needs {self.k} blocks, was given {len(have)}")
        bad = [index for index in have if not 0 <= index < self.k + self.m]
        if bad:
            raise ValueError(f"block indexes {sorted(bad)} are outside 0..{self.k + self.m - 1}")
        indexes = sorted(have)[:self.k]
        blocks = [have[index] for index in indexes]
        length = max(len(block) for block in blocks)
        if indexes == list(range(self.k)):
            return [bytes(block).ljust(length, b"\0") for block in blocks]
        inverse = _invert([self._row(index) for index in indexes])
        return [_combine(inverse[row], blocks, length) for row in range(self.k)]

    def update(self, parity, i, old, new):
        """
        The parity blocks after data block `i` changed from `old` to `new`, without reading any other member.

        Args:
            parity (list of bytes): The `m` parity blocks as they were.
            i (int): The changed member, `0..k-1`.
            old (bytes): Its bytes when the parity was made (empty for a member that was not yet there).
            new (bytes): Its bytes now.

        Returns:
            list of bytes: The new parity, equal to a fresh `encode` up to trailing zero bytes. It is as long as the
                longest of the old parity, `old` and `new`: a longer member extends it, but a shorter one does not
                shrink it (the other members may be as long as the old parity).
        """
        if len(parity) != self.m:
            raise ValueError(f"update takes {self.m} parity blocks, not {len(parity)}")
        if not 0 <= i < self.k:
            raise ValueError(f"member {i} is outside 0..{self.k - 1}")
        length = max(len(old), len(new), *(len(block) for block in parity))
        delta = (int.from_bytes(old, "little") ^ int.from_bytes(new, "little")).to_bytes(length, "little")
        return [(int.from_bytes(block, "little") ^ _scaled(delta, self._matrix[j][i])).to_bytes(length, "little")
                for j, block in enumerate(parity)]
