# tests/test_database_disk.py

"""
The free-space check measures the drive that holds the database's data
directory (asked of the server: `SELECT @@datadir`), never the boot drive
by default, and says "unknown" and why when that drive can't be measured
from this machine.
"""

import os
import shutil
from collections import namedtuple

import pytest

from planetgen.generation import stats as generationStats
from planetgen.generation.stats import DiskSpace, Estimate, Unmeasured, check_disk, database_disk
from planetgen.web.admin_pages import _disk_tile

GB = 1024 ** 3
Usage = namedtuple("Usage", "total used free")


class _Rows:
    def __init__(self, rows):
        self.rows = rows

    def fetchone(self):
        return self.rows[0]

    def fetchall(self):
        return self.rows


class _Conn:
    """A server whose data directory is `datadir` and whose `information_schema.DISKS` is `disks` (or absent)."""

    def __init__(self, datadir, disks=None):
        self.datadir, self.disks = datadir, disks

    def execute(self, sql, *_args):
        if "@@datadir" in sql:
            return _Rows([{"d": self.datadir}])
        if "information_schema.DISKS" in sql:
            if self.disks is None:
                raise RuntimeError("Unknown table 'DISKS'")
            return _Rows(self.disks)
        raise AssertionError(sql)


@pytest.fixture
def two_drives(monkeypatch):
    """The boot drive (`/`) and a second mount at /mnt/big; asking about anything else fails the test."""
    asked = []
    drives = {"/mnt/big": Usage(4000 * GB, 1000 * GB, 3000 * GB), "/": Usage(100 * GB, 90 * GB, 10 * GB)}

    def usage(path):
        asked.append(path)
        return drives["/mnt/big" if str(path).startswith("/mnt/big") else "/"]

    monkeypatch.setattr(shutil, "disk_usage", usage)
    monkeypatch.setattr(os.path, "realpath", lambda path, **_kw: str(path).rstrip("/") or "/")
    monkeypatch.setattr(os.path, "ismount", lambda path: path in ("/", "/mnt/big"))
    return asked


def test_the_data_directory_on_another_drive_is_measured_not_the_boot_drive(two_drives):
    disk = database_disk(_Conn("/mnt/big/mysql/"), "localhost", "planetgen_x")
    assert isinstance(disk, DiskSpace)
    assert (disk.path, disk.mount, disk.free_bytes, disk.total_bytes) == ("/mnt/big/mysql/", "/mnt/big", 3000 * GB, 4000 * GB)
    assert two_drives == ["/mnt/big/mysql/"]
    assert "/mnt/big" in disk.where() and "/mnt/big/mysql/" in disk.where()


def test_a_data_directory_on_the_boot_drive_is_the_boot_drive(two_drives):
    disk = database_disk(_Conn("/var/lib/mysql/"), "127.0.0.1", "planetgen_x")
    assert (disk.mount, disk.free_bytes) == ("/", 10 * GB)


def test_a_symlinked_data_directory_is_resolved_to_the_mount_it_points_at(tmp_path, monkeypatch):
    target = tmp_path / "bigdrive" / "mysql"
    target.mkdir(parents=True)
    link = tmp_path / "var-lib-mysql"
    link.symlink_to(target)
    monkeypatch.setattr(os.path, "ismount", lambda path: str(path) == str(tmp_path / "bigdrive"))
    disk = database_disk(_Conn(str(link)), "localhost", "planetgen_x")
    assert disk.mount == str(tmp_path / "bigdrive")


def test_a_windows_drive_letter_is_the_mount(monkeypatch):
    monkeypatch.setattr(shutil, "disk_usage", lambda path: Usage(2000 * GB, 500 * GB, 1500 * GB))
    monkeypatch.setattr(os.path, "splitdrive", lambda path: ("D:", path[2:]))
    disk = database_disk(_Conn("D:\\MySQL\\Data\\"), "localhost", "planetgen_x")
    assert disk.mount == "D:" + os.sep and disk.free_bytes == 1500 * GB


def test_a_remote_server_is_not_measured_at_a_path_that_only_happens_to_exist_here(two_drives, tmp_path):
    here = tmp_path / "mysql"
    here.mkdir()
    disk = database_disk(_Conn(str(here)), "db.example.com", "planetgen_x")
    assert isinstance(disk, Unmeasured) and "not measurable" in disk.reason
    assert two_drives == []  # the boot drive's numbers were never read


def test_a_remote_server_is_measured_when_its_data_directory_is_really_visible(two_drives, tmp_path, monkeypatch):
    (tmp_path / "planetgen_x").mkdir()
    disk = database_disk(_Conn(str(tmp_path) + "/"), "db.example.com", "planetgen_x")
    assert isinstance(disk, DiskSpace)


def test_a_remote_mariadb_reports_its_own_drive(two_drives):
    disks = [{"p": "/", "t": 100 * GB, "a": 10 * GB}, {"p": "/mnt/big", "t": 4000 * GB, "a": 3000 * GB}]
    disk = database_disk(_Conn("/mnt/big/mysql/", disks), "db.example.com", "planetgen_x")
    assert (disk.source, disk.mount, disk.free_bytes) == ("server report", "/mnt/big", 3000 * GB)
    assert two_drives == []


def test_an_unmeasurable_drive_is_unknown_with_a_reason_and_refuses_nothing():
    result = Estimate(sectors=1, bytes=10 ** 15)
    refusal = check_disk(result, Unmeasured("the database server's data directory (/x) is not measurable here"))
    assert refusal is None and result.disk is None
    assert "not measurable here" in result.disk_note
    assert "Free space not measured" in result.summary() and "nothing is refused" in result.summary()
    assert result.as_dict()["disk_note"] == result.disk_note


def test_the_summary_and_a_refusal_name_the_drive_measured():
    result = Estimate(sectors=1, bytes=2 * GB)
    check_disk(result, DiskSpace("/mnt/big/mysql/", 100 * GB, 6 * GB, mount="/mnt/big"))
    assert "/mnt/big" in result.refusal and "/mnt/big/mysql/" in result.refusal
    assert "on /mnt/big (data directory /mnt/big/mysql/)" in result.summary()


def test_a_server_that_does_not_report_its_data_directory_is_unknown():
    class _Broken:
        def execute(self, *_args):
            raise RuntimeError("no")

    assert isinstance(database_disk(_Broken(), "localhost", "x"), Unmeasured)


def test_the_admin_tile_says_unknown_instead_of_showing_the_boot_drive():
    assert _disk_tile(None) == "unknown"
    assert _disk_tile({"measured": False, "note": "not visible"}).startswith("unknown (not visible")
    shown = _disk_tile({"measured": True, "free_bytes": 3072 * GB, "total_bytes": 4096 * GB, "where": "/mnt/big"})
    assert shown == "3.0 TB of 4.0 TB on /mnt/big"
