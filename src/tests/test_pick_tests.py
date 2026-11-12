"""TEST.722: scripts/pick_tests.py picks test files in priority order, stops at the limit and re-runs only the failed files."""

import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import pick_tests  # noqa: E402


def test_module_names():
    assert pick_tests.module_names("src/planetgen/physics/birth_luminosity.py") == ("planetgen.physics.birth_luminosity", "birth_luminosity")
    assert pick_tests.module_names("src/planetgen/physics/__init__.py") == ("planetgen.physics", "physics")
    assert pick_tests.module_names("scripts/test.py") is None
    assert pick_tests.module_names("src/planetgen/__init__.py") is None


def test_imports_module_counts_the_three_import_forms():
    text = "import planetgen.a.b\nfrom planetgen.a import b, c\nfrom planetgen.a.b import x\nfrom planetgen.a import other\nimport planetgen.a.bc\n"
    assert pick_tests.imports_module(text, "planetgen.a.b") == 3


def test_pick_orders_edited_then_name_then_importers(tmp_path):
    (tmp_path / "test_alpha.py").write_text("import planetgen.m.alpha\n")
    (tmp_path / "test_other.py").write_text("from planetgen.m import alpha\nfrom planetgen.m.alpha import x\nimport planetgen.m.alpha\n")
    (tmp_path / "test_few.py").write_text("from planetgen.m.alpha import x\n")
    (tmp_path / "test_unrelated.py").write_text("import planetgen.m.beta\n")
    changed = ["src/planetgen/m/alpha.py", "src/tests/test_edited.py"]
    picked = pick_tests.pick(changed, tests_dir=tmp_path)
    names = [name.rsplit("/", 1)[-1] for name, _ in picked]
    assert names == ["test_edited.py", "test_alpha.py", "test_few.py", "test_other.py"]
    assert dict(picked)[picked[1][0]].startswith("name holds")


def test_pick_leaves_the_skipped_set_out(tmp_path):
    (tmp_path / "test_alpha.py").write_text("")
    picked = pick_tests.pick(["src/planetgen/m/alpha.py"], tests_dir=tmp_path, skip={str(tmp_path) + "/test_alpha.py"})
    assert picked == []


def test_failed_files_reads_a_pytest_log():
    log = "FAILED src/tests/test_a.py::test_x - assert 1\nERROR src/tests/test_b.py::test_y\nFAILED src/tests/test_a.py::test_z\n3 failed"
    assert pick_tests.failed_files(log) == ["src/tests/test_a.py", "src/tests/test_b.py"]


class _Result:
    def __init__(self, code, out):
        self.returncode, self.stdout = code, out


def test_the_limit_stops_new_batches_and_names_what_is_left():
    files = ["src/tests/test_%d.py" % i for i in range(14)]
    ticks = iter([0, 0, 5])
    calls = []

    def runner(command, **kwargs):
        calls.append(command)
        return _Result(1, "FAILED src/tests/test_0.py::t\n1 failed") if len(calls) == 1 else _Result(0, "6 passed")

    failed, left = pick_tests.run_batches(files, 3, runner=runner, clock=lambda: next(ticks))
    assert failed == ["src/tests/test_0.py"]
    assert len(calls) == 2 and left == files[2 * pick_tests.BATCH:]
