# planetgen/physics/bright_piece_kernel.py

"""
The stretches of a track that clear a luminosity (PERF.232)
===========================================================

`post_main_sequence.Timeline.bright_pieces` finds, for each phase of a star, the stretches of its track (in the share of the
main-sequence lifetime) where it is at least as bright as a threshold. The scatter asks it for every mass of the sampling grid
at every threshold of the luminosity bands, hundreds of thousands of times on a first fill, and the numpy version spends nearly
all of that on the overhead of its two dozen small array operations. This is the same arithmetic, element by element and in the
same order, in one pass over the samples of each phase, so a stretch it returns is the stretch the array version returns, to
the last bit (`test_bright_piece_kernel`).
"""

import numpy as np
from numba import njit


@njit(cache=True)
def bright_ranges(frac, log_l, floor, lows, highs, f_starts, f_ends):
    """
    The stretches of each phase where `log_l` (the track's log luminosity at the samples whose share of the lifetime is `frac`,
    linear between samples) is at least `floor`, joined where one begins no later than 1e-15 after the last ended, as
    `(counts, starts, ends)`: `counts[w]` stretches for phase `w`, one after another in `starts` and `ends`. Phase `w` reads the
    intervals between samples `lows[w]` and `highs[w]` and is cut to `f_starts[w]`..`f_ends[w]`.
    """
    windows = lows.shape[0]
    counts = np.zeros(windows, dtype=np.int64)
    capacity = frac.shape[0] + windows
    starts = np.empty(capacity)
    ends = np.empty(capacity)
    used = 0
    for w in range(windows):
        f_start = f_starts[w]
        f_end = f_ends[w]
        have = False
        piece_start = 0.0
        piece_end = 0.0
        previous_end = 0.0
        for i in range(lows[w], highs[w]):
            e0 = log_l[i] - floor
            e1 = log_l[i + 1] - floor
            above0 = e0 >= 0.0
            above1 = e1 >= 0.0
            if not (above0 or above1):
                continue
            delta = e1 - e0
            if delta == 0.0:
                ratio = -e0 / np.inf
            else:
                ratio = -e0 / delta
            cross = frac[i] + ratio * (frac[i + 1] - frac[i])
            first = frac[i] if above0 else cross
            last = frac[i + 1] if above1 else cross
            start = first if first > f_start else f_start
            end = last if last < f_end else f_end
            if not end > start:
                continue
            if have and not start > previous_end + 1.0e-15:
                piece_end = end
            else:
                if have:
                    starts[used] = piece_start
                    ends[used] = piece_end
                    used += 1
                    counts[w] += 1
                piece_start = start
                piece_end = end
                have = True
            previous_end = end
        if have:
            starts[used] = piece_start
            ends[used] = piece_end
            used += 1
            counts[w] += 1
    return counts, starts[:used], ends[:used]
