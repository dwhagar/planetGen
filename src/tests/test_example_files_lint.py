# tests/test_example_files_lint.py

"""
OPS.31: every example plist, XML and systemd unit file in the repository
parses. A malformed one used to ship unnoticed (OPS.32: a `--` inside the
update plist's comment made it invalid, and the installer's `|| true` hid
it). The files are found with `git ls-files`, so a new example is covered
without touching this test.
"""

import os
import plistlib
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET

import pytest

REPO_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
UNIT_SUFFIXES = (".service", ".timer", ".socket", ".path", ".mount", ".target")


def _tracked(*suffixes, under=("examples/", "scripts/", "docs/")):
    out = subprocess.run(["git", "ls-files"], cwd=REPO_DIR, capture_output=True, text=True, check=False).stdout
    files = [p for p in out.splitlines() if p.endswith(suffixes) and p.startswith(under)]
    if not out:
        pytest.skip("not a git checkout")
    return files


def test_there_are_example_files_to_lint():
    assert _tracked(".plist")
    assert _tracked(*UNIT_SUFFIXES)


@pytest.mark.parametrize("path", _tracked(".plist"))
def test_every_plist_parses(path):
    with open(os.path.join(REPO_DIR, path), "rb") as handle:
        assert isinstance(plistlib.load(handle), dict), path


@pytest.mark.parametrize("path", _tracked(".xml", ".xsd", ".xsl"))
def test_every_xml_file_parses(path):
    ET.parse(os.path.join(REPO_DIR, path))


# systemd-analyze names a command that isn't installed on this machine; the
# examples point at /usr/local/planetgen and /var/lib/planetGen, which
# only exist on a server.
ABSENT_ON_THIS_MACHINE = re.compile(r"Command \S+ is not executable: No such file or directory")


@pytest.mark.skipif(shutil.which("systemd-analyze") is None, reason="needs systemd-analyze")
@pytest.mark.parametrize("path", _tracked(*UNIT_SUFFIXES))
def test_every_systemd_unit_verifies(path, tmp_path):
    unit = tmp_path / os.path.basename(path)
    shutil.copy(os.path.join(REPO_DIR, path), unit)
    proc = subprocess.run(["systemd-analyze", "verify", str(unit)], capture_output=True, text=True, check=False)
    problems = [line for line in (proc.stdout + proc.stderr).splitlines()
                if line.strip() and not ABSENT_ON_THIS_MACHINE.search(line)]
    assert not problems, "\n".join(problems)
