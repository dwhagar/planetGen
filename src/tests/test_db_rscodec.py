# tests/test_db_rscodec.py

"""
DB.43: the Reed-Solomon codec (`planetgen.util.rscodec`) rebuilds any lost blocks exactly, updates parity from one
changed member, and takes members of unequal length through zero padding. No database is needed.
"""

import itertools
import random

import pytest

from planetgen.util import rscodec
from planetgen.util.rscodec import RS

SHAPES = [(16, 1), (16, 2), (32, 2), (8, 3)]


def _blocks(rng, count, length):
    return [bytes(rng.randrange(256) for _ in range(length)) for _ in range(count)]


@pytest.mark.parametrize("k, m", SHAPES)
def test_any_pattern_of_up_to_m_lost_blocks_decodes_exactly(k, m):
    rng = random.Random(4300 + 10 * k + m)
    code = RS(k, m)
    data = _blocks(rng, k, 200)
    group = data + code.encode(data)
    for lost_count in range(0, m + 1):
        for _ in range(25):
            lost = set(rng.sample(range(k + m), lost_count))
            have = {i: block for i, block in enumerate(group) if i not in lost}
            assert code.decode(have) == data, (k, m, sorted(lost))


def test_every_loss_pattern_of_a_small_group_decodes():
    rng = random.Random(7)
    code = RS(4, 3)
    data = _blocks(rng, 4, 33)
    group = data + code.encode(data)
    for lost in itertools.chain.from_iterable(itertools.combinations(range(7), n) for n in range(4)):
        assert code.decode({i: b for i, b in enumerate(group) if i not in lost}) == data


def test_more_blocks_than_needed_are_fine_and_fewer_are_refused():
    rng = random.Random(1)
    code = RS(5, 2)
    data = _blocks(rng, 5, 40)
    group = data + code.encode(data)
    assert code.decode(dict(enumerate(group))) == data
    with pytest.raises(ValueError, match="needs 5"):
        code.decode({i: group[i] for i in (0, 1, 2, 5)})


def test_the_first_parity_block_is_the_plain_xor_of_the_data():
    rng = random.Random(2)
    data = _blocks(rng, 6, 64)
    expected = bytes(a ^ b ^ c ^ d ^ e ^ f for a, b, c, d, e, f in zip(*data))
    assert RS(6, 2).encode(data)[0] == expected
    assert RS(6, 1).encode(data) == [expected]


@pytest.mark.parametrize("k, m", SHAPES)
def test_update_equals_a_fresh_encode(k, m):
    rng = random.Random(900 + k + m)
    code = RS(k, m)
    data = _blocks(rng, k, 150)
    parity = code.encode(data)
    for member in (0, k // 2, k - 1):
        new = _blocks(rng, 1, 150)[0]
        parity = code.update(parity, member, data[member], new)
        data[member] = new
        assert parity == code.encode(data)


def test_a_new_member_in_an_open_group_updates_from_zeros():
    rng = random.Random(3)
    code = RS(8, 2)
    data = [b""] * 8
    parity = [b"", b""]
    for member in range(8):
        data[member] = _blocks(rng, 1, 90 + member)[0]
        parity = code.update(parity, member, b"", data[member])
        fresh = code.encode(data)
        assert [p.rstrip(b"\0") for p in parity] == [p.rstrip(b"\0") for p in fresh]
        assert [len(p) for p in parity] == [len(p) for p in fresh]


def test_unequal_member_lengths_work_through_zero_padding():
    rng = random.Random(5)
    code = RS(6, 2)
    lengths = [10, 300, 0, 299, 1, 57]
    data = [_blocks(rng, 1, n)[0] for n in lengths]
    parity = code.encode(data)
    assert [len(p) for p in parity] == [300, 300]
    group = data + parity
    for lost in ([1, 3], [0, 6], [2, 7], [4]):
        have = {i: b for i, b in enumerate(group) if i not in lost}
        rebuilt = code.decode(have)
        assert [block[:n] for block, n in zip(rebuilt, lengths)] == data
        assert all(len(block) == 300 for block in rebuilt)
        assert all(not block[n:].strip(b"\0") for block, n in zip(rebuilt, lengths))


def test_update_to_a_longer_member_extends_the_parity_and_a_shorter_one_does_not_shrink_it():
    rng = random.Random(6)
    code = RS(4, 2)
    data = _blocks(rng, 4, 20)
    parity = code.encode(data)
    longer = _blocks(rng, 1, 50)[0]
    extended = code.update(parity, 1, data[1], longer)
    data[1] = longer
    assert extended == code.encode(data) and len(extended[0]) == 50
    shorter = longer[:20]
    again = code.update(extended, 1, longer, shorter)
    data[1] = shorter
    assert len(again[0]) == 50
    assert [p.rstrip(b"\0") for p in again] == [p.rstrip(b"\0") for p in code.encode(data)]


def test_a_group_of_the_widest_size_works():
    rng = random.Random(8)
    code = RS(254, 2)
    data = _blocks(rng, 254, 16)
    group = data + code.encode(data)
    have = {i: b for i, b in enumerate(group) if i not in (3, 200)}
    assert code.decode(have) == data


def test_bad_shapes_and_indexes_are_refused():
    for k, m in ((0, 1), (1, 0), (255, 2), (-1, 2)):
        with pytest.raises(ValueError):
            RS(k, m)
    code = RS(3, 1)
    with pytest.raises(ValueError, match="takes 3"):
        code.encode([b"a", b"b"])
    with pytest.raises(ValueError, match="outside"):
        code.decode({0: b"a", 1: b"b", 9: b"c"})
    with pytest.raises(ValueError, match="takes 1"):
        code.update([b"a", b"b"], 0, b"a", b"b")
    with pytest.raises(ValueError, match="member 3"):
        code.update([b"a"], 3, b"a", b"b")


def test_the_field_tables_are_a_group():
    for a in range(1, 256):
        assert rscodec._multiply(a, rscodec._inverse(a)) == 1
        assert rscodec._multiply(a, 1) == a
    assert rscodec._multiply(0x53, 0xCA) == rscodec._multiply(0xCA, 0x53)
    assert rscodec._multiply(2, 0x80) == 0x1D  # x * x^7 reduces by the polynomial
