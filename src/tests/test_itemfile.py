# tests/test_itemfile.py

"""OPS.240: the item file of a windowed stage: what is written is what a window reads back, and a file that is missing, torn or changed fails loudly."""

import pytest

from planetgen.queue import itemfile


def _records(count, suffixes=False):
    return [(shell, shell % 7, -shell if shell % 5 == 0 else shell * 3, 1 + shell % 11, f", {shell * 1.25:.2f} pc from sector_id=9" if suffixes and shell % 3 else "")
            for shell in range(count)]


@pytest.mark.parametrize("suffixes", (False, True))
def test_every_window_reads_back_what_was_written(tmp_path, suffixes):
    records = _records(53, suffixes)
    path = str(tmp_path / "stage.items")
    assert itemfile.write(path, records) == 53
    with itemfile.Reader(path) as reader:
        reader.verify()
        assert reader.count == 53 and reader.has_suffixes is suffixes
        assert reader.window(0, 53) == records
        for start, size in ((0, 7), (7, 7), (13, 1), (52, 5), (53, 3), (20, 0)):
            assert reader.window(start, size) == records[start:start + size]
        stitched = []
        for start in range(0, 53, 8):
            stitched += reader.window(start, 8)
        assert stitched == records


def test_a_stage_with_no_items_is_a_file_of_nothing(tmp_path):
    path = str(tmp_path / "none.items")
    assert itemfile.write(path, []) == 0
    with itemfile.Reader(path) as reader:
        reader.verify()
        assert reader.count == 0 and reader.window(0, 10) == []


def test_a_file_written_again_is_replaced_by_the_same_content(tmp_path):
    path = str(tmp_path / "again.items")
    itemfile.write(path, _records(10, True))
    first = open(path, "rb").read()
    itemfile.write(path, _records(10, True))
    assert open(path, "rb").read() == first
    assert not [name for name in tmp_path.iterdir() if ".part" in name.name]


def test_a_missing_file_fails_and_removing_one_that_is_gone_does_not(tmp_path):
    path = str(tmp_path / "gone.items")
    with pytest.raises(itemfile.ItemFileError, match="missing"):
        itemfile.Reader(path)
    itemfile.remove(path)


@pytest.mark.parametrize("cut", (3, itemfile.HEADER.size + 5, -1))
def test_a_file_cut_short_fails_when_it_is_opened(tmp_path, cut):
    path = tmp_path / "torn.items"
    itemfile.write(str(path), _records(20))
    data = path.read_bytes()
    path.write_bytes(data[:cut])
    with pytest.raises(itemfile.ItemFileError, match="cut short|torn"):
        itemfile.Reader(str(path))


def test_a_changed_byte_fails_the_check(tmp_path):
    path = tmp_path / "flipped.items"
    itemfile.write(str(path), _records(20, True))
    data = bytearray(path.read_bytes())
    data[itemfile.HEADER.size + 30] ^= 0x01
    path.write_bytes(bytes(data))
    with itemfile.Reader(str(path)) as reader, pytest.raises(itemfile.ItemFileError, match="checksum"):
        reader.verify()


def test_a_file_of_another_version_or_kind_fails(tmp_path):
    path = tmp_path / "other.items"
    itemfile.write(str(path), _records(4))
    data = bytearray(path.read_bytes())
    data[4] = 9
    path.write_bytes(bytes(data))
    with pytest.raises(itemfile.ItemFileError, match="not a version"):
        itemfile.Reader(str(path))
    path.write_bytes(b"x" * 64)
    with pytest.raises(itemfile.ItemFileError, match="not a version"):
        itemfile.Reader(str(path))


def test_a_window_past_the_end_fails(tmp_path):
    path = str(tmp_path / "end.items")
    itemfile.write(path, _records(5))
    with itemfile.Reader(path) as reader, pytest.raises(itemfile.ItemFileError, match="not a record 6"):
        reader.window(6, 1)


def test_the_file_of_a_stage_is_in_the_state_folder_under_the_run_and_the_stage(tmp_path, monkeypatch):
    monkeypatch.setenv("PLANETGEN_STATE_DIR", str(tmp_path))
    assert itemfile.path_of("20261118-000000-abcd", "sectors") == str(tmp_path / "stage-items" / "20261118-000000-abcd-sectors.items")
