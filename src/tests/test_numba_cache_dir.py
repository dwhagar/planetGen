"""The compiled kernels cache somewhere writable even when the code tree is read-only to the service user (SEC.37)."""
import os
import subprocess
import sys

import planetgen


def test_an_explicit_numba_cache_dir_is_left_alone(tmp_path):
    assert planetgen.choose_numba_cache_dir({"NUMBA_CACHE_DIR": "/x"}, str(tmp_path)) is None


def test_a_writable_default_folder_is_used(tmp_path, monkeypatch):
    monkeypatch.setattr(planetgen, "NUMBA_CACHE_DEFAULT", str(tmp_path))
    assert planetgen.choose_numba_cache_dir({}, "/nonexistent") == str(tmp_path)


def test_a_read_only_code_tree_gets_a_private_folder(tmp_path, monkeypatch):
    monkeypatch.setattr(planetgen, "NUMBA_CACHE_DEFAULT", str(tmp_path / "missing"))
    monkeypatch.setattr(planetgen, "_writable", lambda folder: folder != "/ro" and folder != "/ro/__pycache__" and os.path.isdir(folder))
    monkeypatch.setattr("tempfile.gettempdir", lambda: str(tmp_path))
    folder = planetgen.choose_numba_cache_dir({}, "/ro")
    assert folder == str(tmp_path / ("planetgen-numba-%d" % os.getuid()))
    assert os.path.isdir(folder)


def test_the_kernel_imports_and_caches_into_the_chosen_folder(tmp_path):
    env = {k: v for k, v in os.environ.items() if k != "NUMBA_CACHE_DIR"}
    env["NUMBA_CACHE_DIR"] = str(tmp_path)
    src = os.path.dirname(os.path.dirname(os.path.abspath(planetgen.__file__)))
    code = f"import sys; sys.path.insert(0, {src!r}); from planetgen.physics import bright_piece_kernel, field_kernel, state_kernel, radiative_kernel; print('ok')"
    out = subprocess.run([sys.executable, "-I", "-c", code], env=env, capture_output=True, text=True)
    assert out.stdout.strip() == "ok", out.stderr
