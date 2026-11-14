"""
Class M tiers, the envelope dial: under a tier the share `abundance_tier.bare_share` of the ecosphere's rocky planets are drawn as late cores (GEN.1323), a bare core with no hydrogen envelope.

Reuses `Planet(late_core=True)` and `planets.apply_envelope`. The decision is one draw in `Planet._draw_body`, stored as `_late_core`, so the settle at the final age neither draws again nor brings the envelope back.
"""

import pytest

from planetgen import tuning
from planetgen.generation import abundance_tier
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import planet_type
from planetgen.util import draw

SEED = 20261114
SYSTEMS = 12
STARS = ("F9V", "G0V", "G2V", "G4V", "G6V", "K0V")


def _systems(tier, count=SYSTEMS, seed=SEED):
    draw.set_run_seed(seed)
    made = []
    for index in range(count):
        config = SystemConfig()
        config.STAR_TYPE = STARS[index % len(STARS)]
        config.PLANETS = True
        if tier is not None:
            config.ABUNDANCE_TIER = tier
        made.append(StarSystem(system_config=config))
    return made


def _fingerprint(systems):
    return [[(getattr(p, "body_type", None), getattr(p, "zone", None), round(float(getattr(p, "mass", 0.0)), 6), round(float(getattr(p, "radius", 0.0) or 0.0), 6), round(float(getattr(p, "envelope_mass_fraction", 0.0) or 0.0), 9), bool(getattr(p, "_late_core", False)))
             for p in system.all_planets] for system in systems]


def _rocky(systems):
    return [p for system in systems for p in system.all_planets if getattr(p, "body_type", None) == "t" and not getattr(p, "is_moon", False)]


def _with_dials(monkeypatch, bare):
    monkeypatch.setitem(tuning.ABUNDANCE_TIER_DIALS, "abundant", (0.0, bare, 0.0))


def test_science_changes_nothing():
    assert _fingerprint(_systems(None)) == _fingerprint(_systems("science"))


def test_a_tier_whose_envelope_dial_is_zero_changes_nothing(monkeypatch):
    plain = _fingerprint(_systems(None))
    _with_dials(monkeypatch, 0.0)
    assert plain == _fingerprint(_systems("abundant"))


def test_at_share_one_no_candidate_keeps_an_envelope_and_the_others_are_not_made_late_cores(monkeypatch):
    _with_dials(monkeypatch, 1.0)
    dialed = _rocky(_systems("abundant", 18))
    candidates = [p for p in dialed if abundance_tier.candidate(p)]
    assert len(candidates) >= 5
    for planet in candidates:
        assert planet._late_core
        assert not getattr(planet, "_captured_envelope", 0.0)
        assert planet_type.has_ground(planet)           # no deep envelope (a trace of hydrogen from the interior is not one)
    # a planet outside the ecosphere is left to the law: some keep an envelope as on science
    others = [p for p in dialed if not abundance_tier.candidate(p) and getattr(p, "_late_core", False) is False]
    assert others


def test_giants_are_never_late_cores(monkeypatch):
    _with_dials(monkeypatch, 1.0)
    for system in _systems("abundant", 12):
        for planet in system.all_planets:
            if getattr(planet, "body_type", None) == "g":
                assert not getattr(planet, "_late_core", False)
