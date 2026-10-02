"""
`generate.py system` never saves a system that misses a forced option
(GEN.49), and refuses contradictory forcing options up front (GEN.50).
"""

import json
import random

import pytest

import generate
from stellarObjects import program_constants
from stellarObjects.systemData import StarSystem
from tests.bughunt_support import run_cli


def test_a_system_that_never_meets_a_forced_option_is_not_saved(tmp_path, monkeypatch, capsys):
    built = []

    class NeverHabitable(StarSystem):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            built.append(self)
            self.unmet_requirements = ["a habitable world"]

    monkeypatch.setattr(generate, "StarSystem", NeverHabitable)
    out = tmp_path / "system.md"
    with pytest.raises(SystemExit) as exc:
        run_cli("system", ["--star-type", "O5V", "+habitable_world", "--output", str(out)])
    assert exc.value.code == 1
    assert len(built) == program_constants.SINGLE_SYSTEM_GENERATION_ATTEMPTS
    assert not out.exists()
    captured = capsys.readouterr()
    assert "nothing was saved" in captured.out + captured.err


def test_a_hot_star_with_a_forced_habitable_world_has_one_or_fails(tmp_path, monkeypatch):
    """O5V failed +habitable_world about 2 systems in 5 and still saved
    them; now every saved system has one, and a run either writes one or
    exits with an error."""
    kept = []
    real = generate.StarSystem

    def build(*args, **kwargs):
        system = real(*args, **kwargs)
        kept.append(system)
        return system

    monkeypatch.setattr(generate, "StarSystem", build)
    for seed in range(8):
        random.seed(seed)
        kept.clear()
        out = tmp_path / f"system-{seed}.md"
        try:
            run_cli("system", ["--star-type", "O5V", "+habitable_world", "--output", str(out)])
        except SystemExit as exc:
            assert exc.code == 1
            assert not out.exists()
            continue
        assert out.exists()
        assert kept[-1].unmet_requirements == []
        assert kept[-1].count_habitable(kept[-1].planets)[0] > 0


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
