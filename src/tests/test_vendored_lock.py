# tests/test_vendored_lock.py

"""OPS.35: the libraries the site ships itself match `static/vendor/VENDORED.json` (scripts/vendor_lock.py)."""

import importlib.util
import json
import os

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _lock_script():
    spec = importlib.util.spec_from_file_location("vendor_lock", os.path.join(REPO, "scripts", "vendor_lock.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_every_shipped_vendor_file_is_in_the_lock_with_its_hash():
    assert _lock_script().differences() == []


def test_the_lock_names_the_version_integrity_licence_and_bundler_of_every_package():
    vendor = _lock_script()
    lock = json.load(open(vendor.LOCK, encoding="utf-8"))
    assert [(p["package"], p["version"]) for p in lock["packages"]] == [(e["package"], e["version"]) for e in vendor.PACKAGES]
    for package in lock["packages"]:
        assert package["integrity"].startswith("sha512-"), package["package"]
        assert package["license"] and package["built_with"] and package["files"], package["package"]


def test_the_notices_name_the_bundler_and_the_licences_bundled_in_shoelace():
    notices = open(os.path.join(_lock_script().VENDOR, "THIRD_PARTY_NOTICES.txt"), encoding="utf-8").read()
    assert _lock_script().ESBUILD in notices
    for word in ("Lit", "BSD-3-Clause", "Floating UI", "VENDORED.json"):
        assert word in notices, word
