"""
GEN.196, PERF.121: `planetgen scatter` draws one kind of scatter. A star pass alone clears and rewrites its own
stars and keeps the other's; the phenomena keep the stars and always have their nucleus; every option belongs to one
kind; the luminosity pass needs the mass pass; each kind is one stage of its own.
"""

import argparse

import pytest

from planetgen.cli import generate as generate_cli
from planetgen.db import store
from planetgen.generation import run_plan
from tests.test_bright_star_scatter import (  # noqa: F401
    EDGE_PC, EXTENTS, E_VALUE, SHAPE, THRESHOLD, SCATTER_KINDS, _clear_args, _scatter_args, light_steps)


def _run(mysql_config, *flags):
    """One `planetgen scatter` run of the flags."""
    args = _scatter_args(mysql_config, *flags)
    run_plan.run_scatter(args)
    return args


def _clear(mysql_config, *flags):
    """One `planetgen clear-scatter` run of the flags."""
    run_plan.run_clear_scatter(_clear_args(mysql_config, *flags))


def _seed(mysql_config):
    store.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_shell_index=8,
                            expected_system_count_at_density_1=E_VALUE, config=mysql_config)
    store.replace_galaxy_shells(EXTENTS, config=mysql_config)


def _star_ids(mysql_config, where):
    conn = store.get_connection(mysql_config)
    try:
        return {row["id"] for row in conn.execute(f"SELECT id FROM bright_stars WHERE {where}").fetchall()}
    finally:
        conn.close()


def _settings(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        return store.bright_star_scatter_settings(conn)[0], store.bright_star_mass_limit(conn)
    finally:
        conn.close()


def _first_scatter(mysql_config, mass="12"):
    """The mass pass then the luminosity pass, as a `galaxy --then-scatter` run lists them."""
    _seed(mysql_config)
    _run(mysql_config, "--mass", "--phenomenon-min-mass", mass)
    _run(mysql_config, "--luminosity")
    return _star_ids(mysql_config, f"initial_mass_sol >= {mass}"), _star_ids(mysql_config, f"initial_mass_sol < {mass}")


def test_clearing_and_drawing_the_luminosity_pass_keeps_the_mass_pass(mysql_config):
    heavy, light = _first_scatter(mysql_config)
    assert heavy and light
    _clear(mysql_config, "--luminosity")
    assert not _star_ids(mysql_config, "initial_mass_sol < 12") and heavy <= _star_ids(mysql_config, "initial_mass_sol >= 12")
    _run(mysql_config, "--luminosity", "--bright-star-min-luminosity", "5000")
    assert heavy <= _star_ids(mysql_config, "initial_mass_sol >= 12")
    assert not light & _star_ids(mysql_config, "initial_mass_sol < 12"), "the lighter stars are new rows"
    assert _settings(mysql_config) == (5000.0, 12.0)


def test_a_scatter_draws_nothing_over_what_is_stored_and_says_what_to_clear(mysql_config):
    heavy, light = _first_scatter(mysql_config)
    for flags, kind in ((("--mass",), "--mass"), (("--luminosity",), "--luminosity"), (("--phenomena",), None),
                        (("--nebula-field",), None)):
        if kind is None:
            continue       # nothing of these is stored yet
        with pytest.raises(RuntimeError, match=f"clear-scatter {kind}"):
            _run(mysql_config, *flags)
    assert _star_ids(mysql_config, "1 = 1") == heavy | light, "the refusals changed nothing"


def test_the_mass_pass_cannot_be_cleared_while_the_bright_stars_are_stored(mysql_config):
    heavy, light = _first_scatter(mysql_config)
    with pytest.raises(RuntimeError, match="clear-scatter --luminosity"):
        _clear(mysql_config, "--mass")
    assert _star_ids(mysql_config, "1 = 1") == heavy | light


def test_a_new_mass_limit_is_clearing_both_passes_and_drawing_both_again(mysql_config):
    _first_scatter(mysql_config)
    _clear(mysql_config, "--luminosity")
    _clear(mysql_config, "--mass")
    assert not _star_ids(mysql_config, "1 = 1")
    _run(mysql_config, "--mass", "--phenomenon-min-mass", "16")
    _run(mysql_config, "--luminosity")
    assert _settings(mysql_config) == (THRESHOLD, 16.0)
    assert _star_ids(mysql_config, "initial_mass_sol >= 16") and _star_ids(mysql_config, "initial_mass_sol < 16")


def test_the_phenomena_leave_the_stars_alone_and_keep_the_nucleus(mysql_config):
    heavy, light = _first_scatter(mysql_config)
    _run(mysql_config, "--phenomena", "--compact-min-mass", "2")
    with pytest.raises(RuntimeError, match="clear-scatter --phenomena"):
        _run(mysql_config, "--phenomena", "--compact-min-mass", "6")
    _clear(mysql_config, "--phenomena")
    _run(mysql_config, "--phenomena", "--compact-min-mass", "6")
    assert _star_ids(mysql_config, "1 = 1") == heavy | light
    conn = store.get_connection(mysql_config)
    try:
        assert store.phenomenon_scatter_settings(conn)[1] == 6.0
        nuclei = conn.execute("SELECT COUNT(*) AS n FROM phenomenon_scatter WHERE kind = 'quasar'"
                              " OR (kind = 'black-hole' AND subtype = 'supermassive')").fetchone()["n"]
    finally:
        conn.close()
    assert nuclei == 1


def test_clearing_a_scatter_with_nothing_stored_says_so_and_changes_nothing(mysql_config, capsys):
    _seed(mysql_config)
    for kind in ("--mass", "--luminosity", "--phenomena", "--nebula-field"):
        _clear(mysql_config, kind)
    assert not _star_ids(mysql_config, "1 = 1")
    assert "has nothing stored" in capsys.readouterr().out


def test_the_luminosity_pass_with_no_mass_pass_stored_says_to_run_the_mass_pass_first(mysql_config):
    _seed(mysql_config)
    with pytest.raises(RuntimeError, match="planetgen scatter --mass"):
        _run(mysql_config, "--luminosity")
    assert not _star_ids(mysql_config, "1 = 1")


def test_the_mass_pass_with_nothing_stored_draws_only_the_mass_pass(mysql_config):
    _seed(mysql_config)
    _run(mysql_config, "--mass", "--phenomenon-min-mass", "12")
    assert _star_ids(mysql_config, "initial_mass_sol >= 12") and not _star_ids(mysql_config, "initial_mass_sol < 12")


def test_a_scatter_names_exactly_one_kind():
    parser = argparse.ArgumentParser(prefix_chars='-+')
    generate_cli.add_scatter_arguments(parser)
    with pytest.raises(SystemExit):
        parser.parse_args([])
    with pytest.raises(SystemExit):
        parser.parse_args(["--mass", "--luminosity"])
    for kind in SCATTER_KINDS:
        assert getattr(parser.parse_args([kind]), kind[2:].replace("-", "_")) is True


def test_every_option_belongs_to_one_kind():
    parser = argparse.ArgumentParser(prefix_chars='-+')
    generate_cli.add_scatter_arguments(parser)
    for flags in (["--luminosity", "--phenomenon-min-mass", "12"], ["--mass", "--compact-min-mass", "2"],
                  ["--phenomena", "--down-to", "10"], ["--mass", "--down-to", "10"], ["--band"],
                  ["--nebula-field", "--phenomenon-min-mass", "12"]):
        with pytest.raises(SystemExit):
            generate_cli.validate_scatter_args(parser.parse_args(flags), parser)
    generate_cli.validate_scatter_args(parser.parse_args(["--band", "--down-to", "10"]), parser)
    generate_cli.validate_scatter_args(parser.parse_args(["--mass", "--phenomenon-min-mass", "12"]), parser)


def test_clear_scatter_names_exactly_one_kind_and_the_band_is_not_one():
    parser = argparse.ArgumentParser(prefix_chars='-+')
    generate_cli.add_clear_scatter_arguments(parser)
    with pytest.raises(SystemExit):
        parser.parse_args([])
    with pytest.raises(SystemExit):
        parser.parse_args(["--mass", "--phenomena"])
    with pytest.raises(SystemExit):
        parser.parse_args(["--band"])


def test_the_plan_takes_none_of_the_scatter_options():
    _parser, parsers = generate_cli.build_parser()
    for flag in ("--phenomenon-min-mass", "--bright-star-min-luminosity", "--compact-min-mass", "--force",
                 "--no-bright-stars", "--bright-stars-only", "--phenomena-only", "--redo-scatters"):
        with pytest.raises(SystemExit):
            parsers["plan"].parse_args([flag, "1"])


def _recording_stats(monkeypatch):
    """The generation stats of the run, kept in memory: `(kind, seconds)` of everything a run records (the starter's
    and the ones the unit that completes a stage records, PERF.106)."""
    from planetgen.generation import run_common
    from planetgen.generation.stats import GenerationStats

    class Recorder(GenerationStats):
        def __init__(self):
            super().__init__()
            self.seen = []

        def record(self, kind, density, seconds, **kwargs):
            self.seen.append((kind, seconds))
            super().record(kind, density, seconds, **kwargs)

        def flush(self):
            pass

    recorder = Recorder()
    monkeypatch.setattr(run_common, "_generation_stats", lambda args: recorder)
    monkeypatch.setattr(run_plan, "_speed_stats", lambda mysql_config: recorder)    # the unit that completes a stage
    return recorder


def test_the_mass_pass_and_the_luminosity_pass_record_their_own_rows(mysql_config, monkeypatch):
    recorder = _recording_stats(monkeypatch)
    _first_scatter(mysql_config)
    kinds = {kind for kind, _seconds in recorder.seen}
    assert {"scatter-mass", "scatter-bright"} <= kinds and "scatter" not in kinds
    assert recorder.pool_rate("scatter-mass", 1) is not None and recorder.pool_rate("scatter-bright", 1) is not None


def test_the_phenomena_shells_record_their_writing_time_apart_from_the_whole(mysql_config, monkeypatch):
    recorder = _recording_stats(monkeypatch)
    _seed(mysql_config)
    run_plan.scatter_phenomena(_scatter_args(mysql_config, "--phenomena"))
    whole = [seconds for kind, seconds in recorder.seen if kind == "phenomena"]
    written = [seconds for kind, seconds in recorder.seen if kind == "phenomena-write"]
    assert whole and len(written) == len(whole)
    assert 0 <= sum(written) <= sum(whole)
