# tests/test_datasets_gate.py

"""
The data gate in front of bulk generation (GEN.590, Boss's "Require the data"): every bulk command and the Sector page's
neighbourhood generation refuse to start, naming the missing datasets and writing nothing, while a downloaded file is
missing; a single system, sector or phenomenon is not gated; the Generate page's second job step is `datasets status`. The
datasets are the offline fixtures of `tests/datasets_support.py`; no test reaches the network.
"""

import os
import subprocess
import sys

import pytest

from planetgen.admin import activity_log
from planetgen.cli import generate as generate_cli
from planetgen.datasets import catalogue, fetch
from planetgen.db import store
from planetgen.generation import run_galaxy
from tests import datasets_support as support

pytestmark = pytest.mark.real_datasets

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


@pytest.fixture(autouse=True)
def cache(tmp_path, monkeypatch):
    """An empty cache folder: nothing has been downloaded."""
    support.reset_caches(tmp_path / "cache", monkeypatch)
    return tmp_path / "cache"


@pytest.fixture
def downloaded(monkeypatch):
    return support.install_all(monkeypatch)


@pytest.fixture
def handlers(monkeypatch):
    """Records which subcommand handlers ran and every activity log line, instead of generating anything."""
    ran, events = [], []
    for name in list(generate_cli._COMMAND_HANDLERS):
        monkeypatch.setitem(generate_cli._COMMAND_HANDLERS, name, lambda args, name=name: ran.append(name))
    monkeypatch.setattr(activity_log, "event", lambda *a, **k: events.append((a, k)))
    return ran, events


def _main(monkeypatch, *argv):
    monkeypatch.setattr(sys, "argv", ["planetgen", *argv])
    generate_cli.main()


def test_nothing_is_missing_once_the_fixtures_are_downloaded(downloaded):
    assert catalogue.missing() == []


def test_the_catalogue_names_each_missing_dataset_with_its_reason():
    missing = catalogue.missing()
    assert [located.dataset.fill for located in missing] == ["SC1", "SC2", "SC3", "AD9", "AD10", "AD11", "SC4", "SC5", "SC6", "SL1", "PS12", "GS1"]
    assert all(located.missing and located.reason for located in missing)


def test_one_altered_file_makes_its_dataset_missing_again(downloaded, cache):
    (cache / "white_dwarf_models" / "sequences.tar").write_bytes(b"not the file")
    fetch._VERIFIED.clear()
    assert [located.dataset.key for located in catalogue.missing()] == ["white_dwarf_models"]


def test_the_gate_names_the_missing_datasets():
    with pytest.raises(run_galaxy.DatasetsMissing) as refused:
        run_galaxy.require_datasets()
    text = str(refused.value)
    assert "nothing was generated" in text
    for key, fill in (("isochrones", "SC1"), ("white_dwarf_models", "SC2"), ("neutron_star_matter", "SC3"), ("ozone_uv", "AD9"), ("hhe_eos", "GS1")):
        assert f"{key} ({fill})" in text


def test_the_gate_passes_with_every_file_there(downloaded):
    run_galaxy.require_datasets()


@pytest.mark.parametrize("argv", [
    ["galaxy", "--shell", "5"],
    ["plan"],
    ["population"],
    ["sector", "--num-sectors", "3"],
])
def test_bulk_runs_refuse_and_write_nothing(handlers, monkeypatch, capsys, argv):
    ran, events = handlers
    with pytest.raises(SystemExit) as exit_info:
        _main(monkeypatch, *argv)
    assert exit_info.value.code == 1
    assert ran == [] and events == []
    out = capsys.readouterr()
    assert "isochrones (SC1)" in out.out + out.err


@pytest.mark.parametrize("argv, command", [
    (["sector"], "sector"),
    (["system", "--output", "-"], "system"),
    (["phenomenon", "--type", "nebula"], "phenomenon"),
    (["check-math"], "check-math"),
])
def test_single_object_runs_are_not_gated(handlers, monkeypatch, argv, command):
    ran, _ = handlers
    _main(monkeypatch, *argv)
    assert ran == [command]


def test_bulk_runs_go_ahead_when_the_data_is_there(downloaded, handlers, monkeypatch):
    ran, _ = handlers
    _main(monkeypatch, "galaxy", "--shell", "5")
    assert ran == ["galaxy"]


def test_neighborhood_generation_refuses_before_touching_the_database(monkeypatch):
    def no_database(*args, **kwargs):
        raise AssertionError("the database was touched")

    monkeypatch.setattr(store, "get_connection", no_database)
    with pytest.raises(run_galaxy.DatasetsMissing, match="isochrones"):
        run_galaxy.generate_sector_neighborhood(1)


def test_neighborhood_estimate_is_not_gated(monkeypatch):
    """An estimate writes nothing, so it isn't refused (it gets as far as the database)."""
    class Reached(Exception):
        pass

    def reached(*args, **kwargs):
        raise Reached

    monkeypatch.setattr(store, "get_connection", reached)
    with pytest.raises(Reached):
        run_galaxy.generate_sector_neighborhood(1, estimate_only=True)


def test_the_neighborhood_route_reports_the_missing_data(monkeypatch):
    """The Sector page's button gets a 409 naming the missing dataset."""
    from planetgen.api import routes
    from planetgen.api.config import Config
    from planetgen.web.app import create_app

    class _Config(Config):
        SECRET_KEY = "test-secret"

    monkeypatch.setattr(routes, "_resolve_requested_write_db_config", lambda: None)
    app = create_app(_Config)
    view = app.view_functions["api.generate_sector_neighborhood_route"]
    while hasattr(view, "__wrapped__"):  # past the admin and rate-limit wrappers
        view = view.__wrapped__
    with app.test_request_context("/api/sectors/1/generate-neighborhood", method="POST", json={}):
        with pytest.raises(routes.ApiError) as error:
            view(1)
    assert error.value.status_code == 409
    assert "isochrones" in str(error.value)


def test_the_status_command_the_generate_page_runs_exits_1_while_data_is_missing(cache):
    env = dict(os.environ, PLANETGEN_DATA_DIR=str(cache))
    result = subprocess.run([sys.executable, "-m", "planetgen.cli.datasets", "status"], capture_output=True, text=True,
                            cwd=REPO, env=env, timeout=120)
    assert result.returncode == 1, result.stdout + result.stderr
    assert "isochrones: missing (SC1)" in result.stdout and "generation will not start" in result.stdout


def test_update_sh_says_that_generation_refuses_without_the_data():
    with open(os.path.join(REPO, "update.sh"), encoding="utf-8") as handle:
        update = handle.read()
    assert "DATASETS_MISSING" in update and "bulk generation refuses to start" in update


def test_the_data_check_answers_when_another_check_holds_the_check_place(monkeypatch, caplog):
    """A read-only data check is never refused for a busy check place: it checks inline and logs the blocker."""
    from planetgen.cli import datasets as datasets_cli
    from planetgen.queue import runrecord, wrapper

    def busy(*args, **kwargs):
        raise runrecord.Busy("run 20261101-000000-abcdef12 is running (Check the database, check kind)")

    monkeypatch.setattr(wrapper, "run", busy)
    monkeypatch.setattr(datasets_cli, "status_lines", lambda chosen, disregard: (["line"], 0))
    monkeypatch.setattr(datasets_cli, "_chosen", lambda keys: [])
    assert datasets_cli.run_status([], False) == 0
    assert "Check the database" in caplog.text
