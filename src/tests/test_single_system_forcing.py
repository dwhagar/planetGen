"""
`planetgen system` never saves a system that misses a forced option
(GEN.49), and refuses contradictory forcing options up front (GEN.50).
"""

import json

import pytest

pass
from planetgen.util import draw
from planetgen import tuning
from planetgen.generation.system import StarSystem
from planetgen.generation import run_system
from tests.bughunt_support import run_cli


def _never_habitable(monkeypatch):
    built = []

    class NeverHabitable(StarSystem):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            built.append(self)
            self.unmet_requirements = ["an asteroid belt"]

    monkeypatch.setattr(run_system, "StarSystem", NeverHabitable)
    return built


def test_a_system_with_no_room_for_a_forced_option_warns_and_keeps_the_last_try(tmp_path, monkeypatch, capsys):
    # GEN.81: the console warns and does what was asked.
    built = _never_habitable(monkeypatch)
    out = tmp_path / "system.md"
    run_cli("system", ["--star-type", "O5V", "+asteroid_belt", "--output", str(out)])
    assert len(built) == tuning.SINGLE_SYSTEM_GENERATION_ATTEMPTS
    assert out.exists()
    captured = capsys.readouterr()
    assert "WARNING:" in captured.out + captured.err
    assert "the last one is kept without it" in captured.out + captured.err


def test_a_system_that_never_meets_a_forced_option_is_not_saved_under_strict(tmp_path, monkeypatch, capsys):
    built = _never_habitable(monkeypatch)
    out = tmp_path / "system.md"
    with pytest.raises(SystemExit) as exc:
        run_cli("system", ["--star-type", "O5V", "+asteroid_belt", "--output", str(out), "--strict"])
    assert exc.value.code == 1
    assert len(built) == tuning.SINGLE_SYSTEM_GENERATION_ATTEMPTS
    assert not out.exists()
    captured = capsys.readouterr()
    assert "No system with an asteroid belt" in captured.out + captured.err


def test_a_hot_star_with_a_habitable_request_is_saved_whether_or_not_it_has_one(tmp_path, monkeypatch):
    """GEN.1130: +habitable_world only raises the odds, so an O5V system that comes out without a habitable world is saved (under `--strict`
    too) and says `unmet`, with no retry of the whole system; one that has a world has it by the request's own conditions."""
    kept = []
    real = StarSystem

    def build(*args, **kwargs):
        system = real(*args, **kwargs)
        kept.append(system)
        return system

    monkeypatch.setattr(run_system, "StarSystem", build)
    for seed in range(4):
        draw.set_run_seed(seed)
        kept.clear()
        out = tmp_path / f"system-{seed}.md"
        run_cli("system", ["--star-type", "O5V", "+habitable_world", "--strict", "--output", str(out)])
        assert out.exists() and len(kept) == 1
        assert kept[-1].unmet_requirements == []
        assert kept[-1].habitable_request in ("aim", "fallback", "unmet")


@pytest.mark.parametrize("spec", [
    {"planets": False, "asteroid_belt": True},
    {"planets": False, "habitable_world": True},
    {"planets": False, "slots": [{"type": "asteroid_belt"}]},
])
def test_contradictory_forcing_in_a_system_file_is_refused(tmp_path, spec, capsys):
    path = tmp_path / "spec.json"
    path.write_text(json.dumps(spec))
    out = tmp_path / "system.md"
    with pytest.raises(SystemExit):
        run_cli("system", ["--system-file", str(path), "--output", str(out)])
    assert not out.exists()


def test_no_planets_and_a_forced_belt_is_refused_on_the_command_line(tmp_path):
    out = tmp_path / "system.md"
    with pytest.raises(SystemExit) as exc:
        run_cli("system", ["-planets", "+asteroid_belt", "--output", str(out)])
    assert exc.value.code == 2
    assert not out.exists()
