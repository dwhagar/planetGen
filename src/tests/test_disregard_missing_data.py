# tests/test_disregard_missing_data.py

"""
GEN.599: a run told in so many words to disregard missing data takes each missing dataset's labelled fallback instead of
refusing (`fetch.set_disregard_missing_data`, `--disregard-missing-data`, the Generate page's box, the API field). Without the
option the refusal stays, and it now names every file to download with its URL and the setup command. The datasets are the
offline fixtures of `tests/datasets_support.py`; no test reaches the network.
"""

import dataclasses
import os
import subprocess
import sys

import pytest

from planetgen.admin import activity_log
from planetgen.cli import generate as generate_cli
from planetgen.datasets import (baraffe_isochrones, catalogue, fetch, harris_clusters, hhe_eos, hitran_cia, ir_ktables, isochrones, neutron_star_matter,
                                noaa_proton_events, ozone_uv, stellar_spectra, white_dwarf_models)
from planetgen.generation import run_galaxy
from planetgen.physics import stellar_abundances as sa
from planetgen.physics import stellar_composition as sc
from planetgen.queue import api_jobs
from planetgen.web import generate_page
from tests import datasets_support as support

pytestmark = pytest.mark.real_datasets

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


@pytest.fixture(autouse=True)
def cache(tmp_path, monkeypatch):
    """An empty cache folder: nothing has been downloaded. The process-wide switch is put back after the test."""
    support.reset_caches(tmp_path / "cache", monkeypatch)
    yield tmp_path / "cache"
    fetch.set_disregard_missing_data(False)
    fetch._STAND_INS.clear()


@pytest.fixture
def mamajek(monkeypatch, cache):
    """The Mamajek sequence downloaded: it has no fallback (GEN.943), so a run told to disregard missing data is still refused while
    it is missing. A test of the gate that lets a run through takes it, so that only datasets with a labelled fallback are missing."""
    installed = support.install_mamajek(monkeypatch)
    monkeypatch.setattr(catalogue, "ALL", tuple(installed if dataset.key == installed.key else dataset for dataset in catalogue.ALL))
    return installed


@pytest.fixture
def handlers(monkeypatch):
    ran = []
    for name in list(generate_cli._COMMAND_HANDLERS):
        monkeypatch.setitem(generate_cli._COMMAND_HANDLERS, name, lambda args, name=name: ran.append(name))
    monkeypatch.setattr(activity_log, "event", lambda *a, **k: None)
    return ran


def _main(monkeypatch, *argv):
    monkeypatch.setattr(sys, "argv", ["planetgen", *argv])
    generate_cli.main()


# --- The refusal (the default) -----------------------------------------------------------------------------------------

def test_a_reader_refuses_by_default_and_names_the_file_its_url_and_the_setup_command():
    with pytest.raises(fetch.DataMissing) as refused:
        sc.star_makeup(sc.NEUTRON_STAR, None, 5.0, 1.4, "ns")
    text = str(refused.value)
    assert neutron_star_matter.COMPOSITIONS in text or "compositions" in text
    assert "https://zenodo.org/records/20085135/files/compositions.zip" in text
    assert "python3 -m planetgen.cli.datasets fetch neutron_star_matter" in text


def test_the_bulk_gate_names_every_file_to_download_with_its_url():
    with pytest.raises(run_galaxy.DatasetsMissing) as refused:
        run_galaxy.require_datasets()
    text = str(refused.value)
    for dataset in (isochrones.DATASET, white_dwarf_models.DATASET, neutron_star_matter.DATASET):
        assert f"python3 -m planetgen.cli.datasets fetch {dataset.key}" in text
        for piece in dataset.pieces:
            assert piece.url in text
    assert "--disregard-missing-data" in text


def test_a_dataset_without_a_fallback_is_refused_whatever_the_run_was_told(monkeypatch):
    bare = dataclasses.replace(isochrones.DATASET, fallback="")
    monkeypatch.setattr(catalogue, "ALL", (bare, white_dwarf_models.DATASET, neutron_star_matter.DATASET))
    with pytest.raises(run_galaxy.DatasetsMissing, match="isochrones"):
        run_galaxy.require_datasets(True)
    fetch.set_disregard_missing_data(True)
    with pytest.raises(fetch.DataMissing):
        fetch.locate(bare).use_stand_in("anything")


# --- The option --------------------------------------------------------------------------------------------------------

def test_the_option_lets_the_gate_through_and_says_what_each_missing_dataset_is_replaced_by(mamajek):
    missing = run_galaxy.require_datasets(True)
    assert [located.dataset.key for located in missing] == ["isochrones", "white_dwarf_models", "neutron_star_matter", "ozone_uv", "ir_ktables", "ir_cia",
                                                        "harris_clusters", "noaa_proton_events", "baraffe_isochrones", "stellar_spectra", "hhe_eos"]
    lines = run_galaxy.disregard_notice(missing)
    assert "--disregard-missing-data" in lines[0]
    assert len(lines) == 11
    for line, dataset in zip(lines[1:], (isochrones.DATASET, white_dwarf_models.DATASET, neutron_star_matter.DATASET, ozone_uv.DATASET,
                                         ir_ktables.DATASET, hitran_cia.DATASET, harris_clusters.DATASET, noaa_proton_events.DATASET,
                                         baraffe_isochrones.DATASET, stellar_spectra.DATASET, hhe_eos.DATASET)):
        assert dataset.key in line and dataset.fallback in line
    assert run_galaxy.disregard_notice([]) == []


def test_a_main_sequence_star_without_the_isochrones_has_no_values_and_says_so():
    with fetch.disregard_missing_data():
        makeup = sc.star_makeup(sc.MAIN_SEQUENCE, 0.0, 4.5, 1.0, "sun")
        assert fetch.stand_ins_used() == [("isochrones", "no main-sequence values")]
    assert makeup.detail["main_sequence"] is None and makeup.detail["main_sequence_lifetime_gy"] is None
    assert any("isochrones are missing" in item for item in makeup.open_items)
    assert makeup.birth is not None, "the composition itself needs no file"
    assert fetch.stand_ins_used() == [] and not fetch.disregarding_missing_data()


def test_a_white_dwarf_without_the_sequences_takes_the_chandrasekhar_radius_and_no_layers():
    with fetch.disregard_missing_data():
        makeup = sc.star_makeup(sc.WHITE_DWARF, -0.2, 8.0, 0.6, "wd", teff_k=20000.0)
    assert makeup.detail["radius_solar"] == pytest.approx(sa.white_dwarf_radius_solar(0.6)) == pytest.approx(0.0127, abs=0.0003)
    assert makeup.detail["layers"] is None
    assert any("white-dwarf sequences are missing" in item for item in makeup.open_items)


def test_a_white_dwarf_outside_the_chandrasekhar_table_has_no_radius_and_does_not_fail():
    with fetch.disregard_missing_data():
        makeup = sc.star_makeup(sc.WHITE_DWARF, -0.2, 8.0, 1.6, "wd", teff_k=20000.0)
    assert makeup.detail["radius_solar"] is None


def test_a_neutron_star_without_the_tables_takes_the_chamel_and_haensel_crust():
    with fetch.disregard_missing_data():
        makeup = sc.star_makeup(sc.NEUTRON_STAR, None, 5.0, 1.4, "ns")
    assert makeup.detail["crust_model"] is None
    assert makeup.detail["outer_crust"][0] == "Fe-56" and makeup.detail["outer_crust"][-1] == "Kr-118"
    assert makeup.detail["neutron_drip_g_cm3"] == sa.NEUTRON_DRIP_DENSITY_G_CM3
    assert any("stand-in sequence" in item for item in makeup.open_items)


def test_the_option_changes_nothing_while_every_file_is_there(monkeypatch):
    support.install_all(monkeypatch)
    cases = [(sc.MAIN_SEQUENCE, 0.0, 4.5, 1.0, "sun", {}), (sc.WHITE_DWARF, -0.2, 8.0, 0.61, "wd", {"teff_k": 60000.0}),
             (sc.NEUTRON_STAR, None, 5.0, 1.4, "ns", {})]
    plain = [sc.star_makeup(*case[:5], **case[5]) for case in cases]
    with fetch.disregard_missing_data():
        again = [sc.star_makeup(*case[:5], **case[5]) for case in cases]
        assert fetch.stand_ins_used() == []
    assert plain == again


def test_a_dataset_missing_one_file_takes_its_fallback_and_the_others_do_not(monkeypatch):
    installed = support.install_all(monkeypatch)
    os.unlink(os.path.join(fetch.folder(installed[2]), neutron_star_matter.COMPOSITIONS))
    fetch._VERIFIED.clear()
    with fetch.disregard_missing_data():
        ns = sc.star_makeup(sc.NEUTRON_STAR, None, 5.0, 1.4, "ns")
        ms = sc.star_makeup(sc.MAIN_SEQUENCE, 0.0, 4.5, 1.0, "sun")
        used = fetch.stand_ins_used()
    assert ns.detail["crust_model"] is None
    assert ms.detail["main_sequence"] is not None
    assert [key for key, _ in used] == ["neutron_star_matter"]


# --- Where the option is given -----------------------------------------------------------------------------------------

@pytest.mark.parametrize("argv", [["galaxy", "--shell", "5"], ["plan"], ["population"]])
def test_a_bulk_command_refuses_without_the_option_and_goes_ahead_with_it(mamajek, handlers, monkeypatch, capsys, argv):
    with pytest.raises(SystemExit) as exit_info:
        _main(monkeypatch, *argv)
    assert exit_info.value.code == 1 and handlers == []
    refusal = capsys.readouterr()
    assert "https://" in refusal.out + refusal.err
    _main(monkeypatch, *argv, "--disregard-missing-data")
    assert handlers == [argv[0]]
    assert fetch.disregarding_missing_data()
    out = capsys.readouterr()
    text = out.out + out.err
    assert "Missing data is disregarded" in text
    for dataset in (isochrones.DATASET, white_dwarf_models.DATASET, neutron_star_matter.DATASET):
        assert dataset.fallback in text


def test_a_command_that_reads_no_data_takes_the_option_and_ignores_it(handlers, monkeypatch, capsys):
    _main(monkeypatch, "check-math", "--disregard-missing-data")
    assert handlers == ["check-math"]
    assert "Missing data is disregarded" not in "".join(capsys.readouterr())


def test_the_option_is_off_in_a_run_that_was_not_told(handlers, monkeypatch):
    fetch.set_disregard_missing_data(True)
    _main(monkeypatch, "system", "--output", "-")
    assert not fetch.disregarding_missing_data()


def test_the_status_command_exits_0_for_a_missing_dataset_with_a_fallback_when_told_to(cache):
    env = dict(os.environ, PLANETGEN_DATA_DIR=str(cache))
    run = lambda *extra: subprocess.run(  # noqa: E731
        [sys.executable, "-m", "planetgen.cli.datasets", "status", *extra], capture_output=True, text=True, cwd=REPO,
        env=env, timeout=120)
    with_fallback = ("isochrones", "white_dwarf_models", "neutron_star_matter")
    refused, told = run(*with_fallback), run(*with_fallback, "--disregard-missing-data")
    assert refused.returncode == 1 and "generation will not start" in refused.stdout
    assert told.returncode == 0, told.stdout + told.stderr
    assert "the run takes a labelled fallback" in told.stdout and "generation will not start" not in told.stdout
    assert "https://mist.science" in told.stdout, "the missing files are still listed with their URLs"
    # the Mamajek sequence has no fallback (GEN.943): told to disregard, the command still exits 1 and says generation will not start
    bare = run("mamajek_dwarf_sequence", "--disregard-missing-data")
    assert bare.returncode == 1 and "generation will not start" in bare.stdout, bare.stdout + bare.stderr


def test_the_neighborhood_job_sets_the_switch_in_its_worker(monkeypatch):
    seen = []
    monkeypatch.setattr(run_galaxy, "generate_sector_neighborhood",
                        lambda sector, radius_ly=None, config=None: seen.append(fetch.disregarding_missing_data()))
    api_jobs.generate_neighborhood(1, 20, None, True)
    api_jobs.generate_neighborhood(1, 20, None)
    assert seen == [True, False]


def test_the_neighborhood_route_takes_the_field(mamajek, monkeypatch):
    from planetgen.api import routes
    from planetgen.api.config import Config
    from planetgen.web.app import create_app

    class _Config(Config):
        SECRET_KEY = "test-secret"

    monkeypatch.setattr(routes, "_resolve_requested_write_db_config", lambda: None)
    asked = []

    class Reached(Exception):
        pass

    def stop(*args, **kwargs):
        asked.append(True)
        raise Reached

    monkeypatch.setattr(run_galaxy, "generate_sector_neighborhood", stop)
    app = create_app(_Config)
    view = app.view_functions["api.generate_sector_neighborhood_route"]
    while hasattr(view, "__wrapped__"):
        view = view.__wrapped__
    with app.test_request_context("/api/sectors/1/generate-neighborhood", method="POST", json={}):
        with pytest.raises(routes.ApiError) as refused:
            view(1)
    assert refused.value.status_code == 409 and not asked
    with app.test_request_context("/api/sectors/1/generate-neighborhood", method="POST",
                                  json={"disregard_missing_data": True}):
        with pytest.raises(Reached):
            view(1)
    assert asked, "the gate let the request through to the estimate"


def test_the_generate_page_hands_the_option_to_every_generating_command(monkeypatch):
    plain = generate_page.build_job("plan", {}, "galaxy")[2]
    ticked = generate_page.build_job("plan", {"disregard_missing_data": "1"}, "galaxy")[2]
    assert all(generate_page.DISREGARD_FLAG not in step["argv"] for step in plain)
    labels = {step["label"]: step["argv"] for step in ticked}
    assert labels[generate_page.DATA_CHECK_LABEL][-2:] == ["status", generate_page.DISREGARD_FLAG]
    generating = [step for step in ticked if step["label"] != generate_page.DATA_CHECK_LABEL]
    assert generating and all(step["argv"][-1] == generate_page.DISREGARD_FLAG for step in generating)


def test_a_reset_or_a_database_check_has_no_data_check_to_tell():
    for action in ("reset", "check_db"):
        kind, _title, steps = generate_page.build_job(action, {"disregard_missing_data": "1", "confirm": "galaxy"}, "galaxy")
        assert all(generate_page.DISREGARD_FLAG not in step["argv"] for step in steps)


def test_the_brown_dwarf_cap_without_the_baraffe_table_is_the_peak_when_told_to_disregard_and_a_refusal_otherwise():
    from planetgen.physics import substellar_ceiling
    with pytest.raises(fetch.DataMissing):
        substellar_ceiling.cap_k(0.05, 1.0)
    with fetch.disregard_missing_data():
        assert substellar_ceiling.cap_k(0.05, 5.0) == substellar_ceiling.HOTTEST_K
        assert ("baraffe_isochrones", "brown dwarf cap at the 2,950 K peak at every age") in fetch.stand_ins_used()
