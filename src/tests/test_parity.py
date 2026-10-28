# tests/test_parity.py

"""
DB.46: the parity file (`db/parity_file.py`), the groups and their updates (`db/parity.py`), the sector export
(`db/sector_export.py`) and the hooks that keep the parity current as sectors are saved, edited and deleted.
"""

import os
import random
import struct
import types

import pytest

from planetgen.db import edits, parity, parity_file as pfile, sector_export, storage_checksum as sc, store
from planetgen.util import rscodec

from tests.test_storage_checksum import _sector_with_system

K, M = 4, 2


def _blocks(rng, count, low=900, high=1200):
    return [rng.randbytes(rng.randint(low, high)) for _ in range(count)]


def _fresh(blocks, k=K, m=M):
    """The parity a fresh encode of `blocks` (zero-padded to `k`) gives, as a group record holds it."""
    padded = list(blocks) + [b""] * (k - len(blocks))
    longest = max(len(block) for block in padded)
    return [block[:longest] for block in rscodec.RS(k, m).encode(padded)]


@pytest.fixture
def pstore(tmp_path):
    made = parity.ParityStore(str(tmp_path / "t.parity"), K, M)
    made.file.sync = False
    return made


def _group_of(pstore, sector_id):
    with pstore.file.locked():
        group_id, slot, _held = pstore.member(sector_id)
        return pstore.file.read(group_id), slot


# -- the file ---------------------------------------------------------------------------------------------------

def test_a_record_survives_the_round_trip_and_its_crc_catches_a_flipped_byte():
    rng = random.Random(1)
    record = pfile.Record(7, 3, K, M, [pfile.Member(5, 10, b"\1" * 32)] + [pfile.Member()] * (K - 1),
                          [rng.randbytes(10) for _ in range(M)])
    data = record.encode()
    again = pfile.decode(data)
    assert (again.group_id, again.seq, again.k, again.m) == (7, 3, K, M)
    assert again.members == record.members and again.parity == record.parity
    flipped = bytearray(data)
    flipped[-9] ^= 1
    with pytest.raises(pfile.RecordDamaged):
        pfile.decode(bytes(flipped))


def test_the_last_record_of_a_group_wins_and_a_scan_of_the_file_finds_the_same_groups(pstore, tmp_path):
    rng = random.Random(2)
    first = _blocks(rng, 3)
    for sector_id, block in zip((10, 20, 30), first):
        pstore.put(sector_id, block, bytes([sector_id]) * 32)
    pstore.put(20, b"new twenty", b"\x99" * 32, read_member=lambda other: (first[0 if other == 10 else 2], bytes([other]) * 32))
    again = pfile.ParityFile(pstore.file.path)         # a second process, scanning
    with again.locked():
        assert again.sector_group == pstore.file.sector_group
        assert set(again.entries) == set(pstore.file.entries)
        record = again.read(next(iter(again.entries)))
    assert [m.sector_id for m in record.members if m.sector_id] == [10, 20, 30]
    assert [m.checksum[0] for m in record.members[:3]] == [10, 0x99, 30]


def test_a_torn_tail_loses_nothing_and_is_written_over(pstore):
    rng = random.Random(3)
    blocks = _blocks(rng, 3)
    for sector_id, block in zip((100, 200, 300), blocks):
        pstore.put(sector_id, block, bytes([sector_id % 256]) * 32)
    path = pstore.file.path
    with pstore.file.locked():
        good = pstore.file.valid_end
    with open(path, "ab") as handle:       # a record cut off half way
        handle.write(pfile.Record(9, 1, K, M, [pfile.Member(1, 5, b"\0" * 32)] * K, [b"abcde"] * M).encode()[:60])
    again = pfile.ParityFile(path)
    with again.locked():
        assert again.valid_end == good and again.damaged_at is None
        assert set(again.sector_group) == {100, 200, 300}
    again.sync = False
    with again.locked():
        record = again.read(next(iter(again.entries)))
        record.seq += 1
        again.append(record)                   # overwrites the torn bytes
        assert os.path.getsize(path) == again.valid_end
    third = pfile.ParityFile(path)
    with third.locked():
        assert set(third.sector_group) == {100, 200, 300}


def test_a_torn_last_record_with_a_bad_crc_is_ignored(pstore):
    rng = random.Random(4)
    pstore.put(1, rng.randbytes(500), b"\1" * 32)
    pstore.put(2000, rng.randbytes(500), b"\2" * 32)
    with pstore.file.locked():
        end = pstore.file.valid_end
    with open(pstore.file.path, "r+b") as handle:
        handle.seek(end - 3)
        handle.write(b"\xff\xff\xff")
    again = pfile.ParityFile(pstore.file.path)
    with again.locked():
        assert again.damaged_at is None
        assert again.valid_end < end            # the second write is dropped, the first stands


def test_damage_behind_the_last_record_is_reported_and_appends_are_refused_until_it_is_compacted(pstore):
    rng = random.Random(5)
    for sector_id in (1, 2000, 4000, 6000):
        pstore.put(sector_id, rng.randbytes(500), bytes([sector_id % 256]) * 32)
    with open(pstore.file.path, "r+b") as handle:       # a head damaged in the middle of the log
        handle.seek(0)
        handle.write(b"XXXX")
    again = pfile.ParityFile(pstore.file.path)
    again.sync = False
    with again.locked():
        assert again.damaged_at == 0
        with pytest.raises(pfile.ParityFileError):
            again.append(pfile.Record(1, 1, K, M, [pfile.Member()] * K, [b""] * M))


def test_compaction_keeps_only_the_latest_record_of_each_group(pstore):
    rng = random.Random(6)
    blocks = _blocks(rng, 4)
    for sector_id, block in zip((1, 2, 3, 4), blocks):
        pstore.put(sector_id, block, bytes([sector_id]) * 32)
    with pstore.file.locked():
        before = os.path.getsize(pstore.file.path)
        live = pstore.file.live_bytes()
        assert before > live
        assert pstore.file.compact() == (1, 0)
    assert os.path.getsize(pstore.file.path) == live
    again = pfile.ParityFile(pstore.file.path)
    with again.locked():
        assert set(again.sector_group) == {1, 2, 3, 4}
        assert again.read(next(iter(again.entries))).parity == _fresh(blocks)[:M]


# -- the groups -------------------------------------------------------------------------------------------------

def test_a_joining_member_gives_the_parity_a_fresh_encode_would(pstore):
    rng = random.Random(7)
    blocks = _blocks(rng, 3)
    for sector_id, block in zip((10, 20, 30), blocks):
        assert pstore.put(sector_id, block, bytes([sector_id]) * 32) == "joined"
    record, _slot = _group_of(pstore, 10)
    assert record.parity == _fresh(blocks)
    assert [m.length for m in record.members[:3]] == [len(b) for b in blocks]


@pytest.mark.parametrize("through", ["delta", "recomputed"])
def test_a_changed_member_keeps_the_parity_equal_to_a_fresh_encode(pstore, through):
    rng = random.Random(8)
    blocks = _blocks(rng, 4)
    sums = {10: b"\1" * 32, 20: b"\2" * 32, 30: b"\3" * 32, 40: b"\4" * 32}
    for sector_id, block in zip((10, 20, 30, 40), blocks):
        pstore.put(sector_id, block, sums[sector_id])
    longer = rng.randbytes(len(blocks[1]) + 400)
    if through == "delta":
        mode = pstore.put(20, longer, b"\x22" * 32, before=(blocks[1], sums[20]))
    else:
        others = {10: blocks[0], 30: blocks[2], 40: blocks[3]}
        mode = pstore.put(20, longer, b"\x22" * 32, read_member=lambda other: (others[other], sums[other]))
    assert mode == through
    blocks[1] = longer
    record, slot = _group_of(pstore, 20)
    assert record.parity == _fresh(blocks)
    assert record.members[slot].checksum == b"\x22" * 32 and record.members[slot].length == len(longer)


def test_a_before_image_that_is_not_the_one_in_the_record_is_not_used_for_a_delta(pstore):
    rng = random.Random(9)
    blocks = _blocks(rng, 2)
    pstore.put(10, blocks[0], b"\1" * 32)
    pstore.put(20, blocks[1], b"\2" * 32)
    with pytest.raises(parity.NeedBefore):
        pstore.put(20, b"new", b"\3" * 32, before=(blocks[1], b"\xee" * 32))


def test_a_sector_taken_out_leaves_the_parity_of_the_others(pstore):
    rng = random.Random(10)
    blocks = _blocks(rng, 4)
    sums = [bytes([n]) * 32 for n in (1, 2, 3, 4)]
    for sector_id, block, checksum in zip((10, 20, 30, 40), blocks, sums):
        pstore.put(sector_id, block, checksum)
    assert pstore.release(20, before=(blocks[1], sums[1])) == "delta"
    record, _slot = _group_of(pstore, 10)
    assert record.parity == _fresh([blocks[0], b"", blocks[2], blocks[3]])
    assert pstore.release(30, read_member=lambda other: ({10: blocks[0], 40: blocks[3]}[other], sums[other // 10 - 1])) == "recomputed"
    record, _slot = _group_of(pstore, 10)
    assert record.parity == _fresh([blocks[0], b"", b"", blocks[3]])
    assert pstore.release(999) == "absent"
    assert 20 not in pstore.members() and 30 not in pstore.members()


def test_a_free_slot_is_taken_by_the_next_sector_without_reading_any_other(pstore):
    rng = random.Random(11)
    blocks = _blocks(rng, 3)
    for sector_id, block in zip((10, 20, 30), blocks):
        pstore.put(sector_id, block, bytes([sector_id]) * 32)
    pstore.release(20, before=(blocks[1], bytes([20]) * 32))
    newcomer = rng.randbytes(1000)
    assert pstore.put(70, newcomer, b"\x07" * 32) == "joined"
    record, slot = _group_of(pstore, 70)
    assert slot == 1 and record.parity == _fresh([blocks[0], newcomer, blocks[2]])


def test_a_lost_member_is_rebuilt_from_the_others_and_the_parity(pstore):
    rng = random.Random(12)
    blocks = _blocks(rng, 4)
    for sector_id, block in zip((10, 20, 30, 40), blocks):
        pstore.put(sector_id, block, bytes([sector_id]) * 32)
    record, _slot = _group_of(pstore, 10)
    code = rscodec.RS(record.k, record.m)
    have = {0: blocks[0], 1: blocks[1], 3: blocks[3], K + 1: record.parity[1], K: record.parity[0]}
    have.pop(1)
    have.pop(3)                                         # two lost: the code's whole strength
    rebuilt = code.decode({**have, 0: blocks[0], 2: blocks[2]})
    assert rebuilt[1][:len(blocks[1])] == blocks[1] and rebuilt[3][:len(blocks[3])] == blocks[3]


def test_members_far_apart_in_sector_id_share_a_group_and_near_ones_do_not():
    index = parity.GroupIndex()
    members = lambda ids, k=4: tuple(pfile.Member(i, 1000, b"\0" * 32) for i in ids) + (pfile.Member(),) * (k - len(ids))
    entry = lambda ids: pfile.Entry(0, 0, 1, 4, 2, 1000, members(ids))
    index.changed(1, None, entry([100]))
    index.changed(2, None, entry([700]))
    # a big galaxy (window 256): 150 is within 256 of 100, 400 is not within 256 of 100 or 700
    assert index.choose(150, 1000, 100_000) == 2
    assert index.choose(400, 1000, 100_000) == 1            # the fullest first; both have one member, the lowest id wins
    assert index.choose(450, 1000, 100_000) == 1
    index.changed(3, None, entry([100, 700, 1300]))
    assert index.choose(2000, 1000, 100_000) == 3           # the fullest open group
    assert index.choose(120, 1000, 100_000) is None or index.choose(120, 1000, 100_000) not in (1, 3)
    assert index.choose(400, 2_000_000, 100_000) is None    # another size class
    assert parity.window(1000) == 10 and parity.window(10_000_000) == 256 and parity.window(5) == 1


def test_the_size_class_grows_by_a_quarter():
    assert parity.size_class(1) == parity.size_class(2000) == 0
    assert parity.size_class(2001) == 1 and parity.size_class(2500) == 1 and parity.size_class(2501) == 2
    assert parity.size_class(200_000) == 21


def test_the_overhead_on_a_thousand_sectors_of_one_size_is_about_seven_percent(tmp_path):
    rng = random.Random(13)
    made = parity.ParityStore(str(tmp_path / "big.parity"), 32, 2)
    made.file.sync = False
    data = 0
    for sector_id in range(1, 1001):
        block = rng.randbytes(rng.randint(28_000, 32_000))
        made.put(sector_id * 3, block, os.urandom(32))
        data += len(block)
    with made.file.locked():
        groups = len(made.file.entries)
        live = sum(made.file.read(group_id).parity_length * 2 for group_id in made.file.entries)
        assert made.file.compact()[1] == 0
    assert groups < 40
    assert live / data < 0.09
    assert os.path.getsize(made.file.path) < data * 0.11


# -- the decision table ----------------------------------------------------------------------------------------

@pytest.mark.parametrize("rows_match, record, stored, ok, expected", [
    (True, b"a", b"a", True, parity.HEALTHY),
    (True, b"old", b"new", True, parity.RECORD_STALE),
    (False, b"a", b"a", True, parity.REBUILDABLE),
    (False, b"old", b"new", True, parity.LOST),
    (True, b"a", b"a", False, parity.PARITY_DAMAGED),
    (True, None, b"a", True, parity.NOT_PROTECTED),
])
def test_the_stale_record_table(rows_match, record, stored, ok, expected):
    assert parity.classify(rows_match, record, stored, record_ok=ok) == expected


# -- the export -------------------------------------------------------------------------------------------------

@pytest.fixture
def conn(mysql_config):
    connection = store.get_connection(mysql_config)
    yield connection
    connection.close()


def test_an_export_reads_back_into_the_rows_it_was_made_from(conn):
    sector_id, _system_id = _sector_with_system(conn)
    tables = sc.sector_rows(conn, sector_id)
    data, checksum = sector_export.export_sector(conn, sector_id)
    again = sector_export.parse(data)
    assert checksum == sc.checksum(conn, sector_id) == sc.checksum_of(again)
    assert [name for name, _c, _r in again] == [name for name, _c, _r in tables]
    for (name, columns, rows), (_n, again_columns, again_rows) in zip(tables, again):
        assert columns == again_columns
        assert [tuple(row[c] for c in columns) for row in rows] == [tuple(row[c] for c in columns) for row in again_rows], name
    assert sector_export.export_sector(conn, sector_id)[0] == data        # the same bytes every time
    with pytest.raises(sector_export.ExportError):
        sector_export.parse(data[:-5])
    with pytest.raises(sector_export.ExportError):
        sector_export.parse(b"nope")


# -- the database side ------------------------------------------------------------------------------------------

@pytest.fixture
def on(tmp_path, monkeypatch):
    """Parity switched on, the file in `tmp_path`, groups of four sectors with two parity blocks."""
    settings = types.SimpleNamespace(enabled=True, dir=str(tmp_path), group_size=K, parity_blocks=M)
    monkeypatch.setattr(parity, "get_settings", lambda: types.SimpleNamespace(parity=settings))
    monkeypatch.setattr(parity, "size_class", lambda length: 0)         # real sectors differ in size; one class puts them together
    parity._STORES.clear()
    made = parity.store()
    made.file.sync = False
    yield made
    parity._STORES.clear()


def _group_exports(conn, ids):
    return [sector_export.export_sector(conn, sector_id)[0] for sector_id in ids]


def test_a_saved_sector_is_protected_when_the_save_commits(conn, on):
    sector_id, _system_id = _sector_with_system(conn)
    row = conn.execute("SELECT parity_dirty AS d, storage_checksum AS c FROM sectors WHERE id = ?", (sector_id,)).fetchone()
    assert row["d"] == 0 and bytes(row["c"]) == sc.checksum(conn, sector_id)
    with on.file.locked():
        group_id, _slot, held = on.member(sector_id)
        assert held.checksum == bytes(row["c"])
        assert on.file.read(group_id).parity == _fresh(_group_exports(conn, [sector_id]))
    assert parity.member_state(conn, sector_id, on) == parity.HEALTHY
    assert sc.check_sector(conn, sector_id) == sc.OK


def test_an_edit_keeps_the_parity_equal_to_a_fresh_encode(conn, on):
    first, _ = _sector_with_system(conn, "First")
    second, system_id = _sector_with_system(conn, "Second")
    with on.file.locked():
        assert on.member(first)[0] == on.member(second)[0]          # one size class, a small galaxy: one group
    with conn:
        conn.execute("UPDATE star_systems SET name = ? WHERE id = ?", ("Edited Name", system_id))
        sc.mark_system(conn, system_id)
    assert sc.check_sector(conn, second) == sc.OK
    with on.file.locked():
        group_id = on.member(second)[0]
        record = on.file.read(group_id)
    members = [m.sector_id for m in record.members if m.sector_id]
    assert sorted(members) == sorted([first, second])
    exports = {sector_id: sector_export.export_sector(conn, sector_id) for sector_id in members}
    assert record.parity == _fresh([exports[m.sector_id][0] if m.sector_id else b"" for m in record.members], K, M)
    assert all(m.checksum == exports[m.sector_id][1] for m in record.members if m.sector_id)


def test_a_deleted_sector_leaves_the_parity_of_the_rest(conn, on):
    first, _ = _sector_with_system(conn, "First")
    second, _ = _sector_with_system(conn, "Second")
    with conn:
        assert edits.delete_sector_with_contents(conn, first) is not None
    assert first not in on.members() and second in on.members()
    with on.file.locked():
        record = on.file.read(on.member(second)[0])
    slots = [m.sector_id for m in record.members]
    assert record.parity == _fresh([sector_export.export_sector(conn, i)[0] if i else b"" for i in slots], K, M)


def test_a_damaged_member_is_not_recomputed_over(conn, on):
    first, _ = _sector_with_system(conn, "First")
    second, system_id = _sector_with_system(conn, "Second")
    planet = conn.execute("SELECT p.id AS id, p.mass_kg AS mass FROM planets p JOIN star_systems s ON s.id = p.star_system_id "
                          "WHERE s.sector_id = ? ORDER BY p.id LIMIT 1", (first,)).fetchone()
    conn.execute("UPDATE planets SET mass_kg = ? WHERE id = ?", (planet["mass"] * (1 + 2 ** -52), planet["id"]))
    conn.commit()
    assert sc.check_sector(conn, first) == sc.DAMAGED
    with on.file.locked():
        before = on.file.valid_end
    with conn:
        conn.execute("UPDATE star_systems SET name = ? WHERE id = ?", ("Edited", system_id))
        sc.mark_system(conn, system_id)
    assert conn.execute("SELECT parity_dirty AS d FROM sectors WHERE id = ?", (second,)).fetchone()["d"] == 1
    with on.file.locked():
        assert on.file.valid_end == before                          # nothing was written over the damage
    assert parity.member_state(conn, first, on) == parity.REBUILDABLE


def test_a_record_behind_the_database_is_stale_and_is_made_again_not_called_damage(conn, on):
    first, _ = _sector_with_system(conn, "First")
    second, system_id = _sector_with_system(conn, "Second")
    # A crash after the database commit and before the parity record: the rows and the checksum are new, the record old.
    conn.execute("UPDATE star_systems SET name = ? WHERE id = ?", ("Changed", system_id))
    conn.commit()
    with conn:
        sc.seal(conn, second)
    assert parity.member_state(conn, second, on) == parity.RECORD_STALE
    assert sc.check_sector(conn, second) == sc.OK                  # not damage as far as the rows go
    with conn:
        conn.execute("UPDATE sectors SET parity_dirty = 1 WHERE id = ?", (second,))
    assert parity.protect_sector(conn, second, on) == parity.PROTECTED
    assert parity.member_state(conn, second, on) == parity.HEALTHY


def test_a_sector_whose_record_already_holds_its_rows_only_has_its_flag_cleared(conn, on):
    sector_id, _ = _sector_with_system(conn)
    with conn:
        conn.execute("UPDATE sectors SET parity_dirty = 1 WHERE id = ?", (sector_id,))
    assert parity.protect_sector(conn, sector_id, on) == parity.CURRENT
    assert sc.check_sector(conn, sector_id) == sc.OK


def test_a_sector_rebuilt_from_its_group_has_the_rows_and_the_checksum_it_had(conn, on):
    ids = [_sector_with_system(conn, f"S{n}")[0] for n in range(3)]
    exports = {i: sector_export.export_sector(conn, i) for i in ids}
    with on.file.locked():
        record = on.file.read(on.member(ids[0])[0])
    code = rscodec.RS(record.k, record.m)
    lost = ids[1]
    have = {slot: exports[m.sector_id][0] for slot, m in enumerate(record.members) if m.sector_id and m.sector_id != lost}
    have.update({slot: b"" for slot, m in enumerate(record.members) if not m.sector_id})
    have[record.k] = record.parity[0]
    slot = [m.sector_id for m in record.members].index(lost)
    rebuilt = code.decode({**have, record.k: record.parity[0]})[slot][:record.members[slot].length]
    assert rebuilt == exports[lost][0]
    assert sc.checksum_of(sector_export.parse(rebuilt)) == exports[lost][1] == record.members[slot].checksum


def test_refresh_makes_every_sector_current_and_lets_go_of_the_ones_that_are_gone(conn, on, monkeypatch):
    ids = [_sector_with_system(conn, f"R{n}")[0] for n in range(3)]
    with on.file.locked():
        assert set(ids) <= set(on.file.sector_group)
    monkeypatch.setattr(parity, "get_settings", lambda: types.SimpleNamespace(parity=types.SimpleNamespace(
        enabled=False, dir=on.file.path.rsplit("/", 1)[0], group_size=K, parity_blocks=M)))
    parity._STORES.clear()
    on = parity.store()
    on.file.sync = False
    extra, _ = _sector_with_system(conn, "Late")               # parity off: unsealed and not in the file
    with conn:
        conn.execute("DELETE FROM star_systems WHERE sector_id = ?", (ids[2],))
        conn.execute("DELETE FROM sectors WHERE id = ?", (ids[2],))     # a sector gone behind the file's back
    said = []
    counts = parity.refresh(conn, say=said.append)
    assert counts == {parity.PROTECTED: 1, parity.RELEASED: 1}
    assert extra in on.members() and ids[2] not in on.members()
    assert parity.refresh(conn, say=said.append) == {}


def test_nothing_is_written_while_parity_is_off(conn, tmp_path, monkeypatch):
    settings = types.SimpleNamespace(enabled=False, dir=str(tmp_path / "off"), group_size=K, parity_blocks=M)
    monkeypatch.setattr(parity, "get_settings", lambda: types.SimpleNamespace(parity=settings))
    sector_id, _ = _sector_with_system(conn)
    assert not os.path.exists(tmp_path / "off")
    assert sc.check_sector(conn, sector_id) == sc.UNPROTECTED


def test_the_settings_default_to_off_with_two_parity_blocks_in_groups_of_32():
    from planetgen.util.settings import Settings
    defaults = Settings().parity
    assert (defaults.enabled, defaults.group_size, defaults.parity_blocks, defaults.dir) == (False, 32, 2, "")


def test_struct_sizes_match_the_documented_layout():
    assert pfile.HEADER_SIZE == 32 and pfile.SLOT_SIZE == 44 and struct.calcsize("<I") == 4


def test_the_refresh_command_is_registered_and_makes_the_parity_current(conn, on, mysql_config, monkeypatch, capsys):
    from planetgen.cli import generate
    args = generate.build_parser()[0].parse_args(["refresh-parity", "--compact", "--limit", "5"])
    assert args.command == "refresh-parity" and args.compact and args.limit == 5
    assert generate._COMMAND_HANDLERS["refresh-parity"] is generate.run_refresh_parity
    assert "refresh-parity" in generate.READ_ONLY_COMMANDS and generate.DATA_CHANGES["refresh-parity"] == ()
    sector_id, _ = _sector_with_system(conn)
    with on.file.locked():
        before = on.file.valid_end
    with conn:
        conn.execute("UPDATE sectors SET parity_dirty = 1 WHERE id = ?", (sector_id,))
    monkeypatch.setattr(store, "mysql_config_from_args", lambda parsed: mysql_config)
    generate.run_refresh_parity(args)
    assert sc.check_sector(conn, sector_id) == sc.OK
    with on.file.locked():
        assert on.file.valid_end == before                 # the record already held these rows: only the flag was cleared
        assert os.path.getsize(on.file.path) == on.file.valid_end
    assert "compacted to 1 group record" in capsys.readouterr().out
