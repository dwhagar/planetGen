# tests/test_draw_bulk.py

"""
PERF.91: a stream's bulk draws (`Stream.random_array`, `Stream.randbits_list`) give the same numbers one call at a
time would and leave the stream in the same state, so no galaxy changes and nothing is reseeded.
"""

import pytest

from planetgen.util import draw

COUNTS = [0, 1, 7, draw.BULK_MIN - 1, draw.BULK_MIN, draw.BULK_MIN + 1, 5000]


@pytest.mark.parametrize("count", COUNTS)
@pytest.mark.parametrize("seed", ["a stream", 12345])
def test_a_bulk_draw_is_the_same_floats_and_leaves_the_same_state(count, seed):
    bulk, loop = draw.Stream(seed), draw.Stream(seed)
    assert bulk.random_array(count).tolist() == [loop.random() for _ in range(count)]
    assert bulk.getstate() == loop.getstate()
    assert bulk.random() == loop.random()


@pytest.mark.parametrize("bits", [1, 31, 32, 33, 63, 64])
@pytest.mark.parametrize("count", [3, draw.BULK_MIN // 2, draw.BULK_MIN, 4000])
def test_bulk_random_bits_are_the_same_values_and_leave_the_same_state(bits, count):
    bulk, loop = draw.Stream("bits"), draw.Stream("bits")
    values = bulk.randbits_list(bits, count)
    assert values == [loop.getrandbits(bits) for _ in range(count)]
    assert all(isinstance(value, int) and 0 <= value < 1 << bits for value in values)
    assert bulk.getstate() == loop.getstate()


def test_bits_beyond_64_are_drawn_the_ordinary_way():
    bulk, loop = draw.Stream(9), draw.Stream(9)
    assert bulk.randbits_list(100, 1200) == [loop.getrandbits(100) for _ in range(1200)]
    assert bulk.getstate() == loop.getstate()


def test_bulk_draws_between_ordinary_ones_keep_the_sequence():
    bulk, loop = draw.Stream("mixed"), draw.Stream("mixed")
    seen_bulk = [bulk.random() for _ in range(5)]
    seen_loop = [loop.random() for _ in range(5)]
    seen_bulk += bulk.random_array(3000).tolist() + bulk.randbits_list(63, 2000) + [bulk.gauss() for _ in range(3)]
    seen_loop += [loop.random() for _ in range(3000)] + [loop.getrandbits(63) for _ in range(2000)] \
        + [loop.gauss() for _ in range(3)]
    seen_bulk += [bulk.random() for _ in range(5)]
    seen_loop += [loop.random() for _ in range(5)]
    assert seen_bulk == seen_loop and bulk.getstate() == loop.getstate()


def test_a_bulk_draw_is_faster_than_the_loop_above_the_threshold():
    """Not a timing gate: the bulk path is the one taken from BULK_MIN draws on."""
    import time
    count = 200_000
    started = time.perf_counter()
    draw.Stream("t").random_array(count)
    bulk = time.perf_counter() - started
    stream = draw.Stream("t")
    started = time.perf_counter()
    [stream.random() for _ in range(count)]
    assert bulk < (time.perf_counter() - started)
