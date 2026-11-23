# tests/test_class_m_tier_desert.py

"""
GEN.1721, second slice: a desert world's small sea (a halite-saturated sea, GEN.1730) is a living sea for the oxygen its life makes and not a brine seep, and under a Class M abundance tier a share of the candidates has a
quiet interior whose volcanic gases take a third of the oxygen. `science` changes nothing: no draw, no field on the body, the same source and sink.
"""

from types import SimpleNamespace

import pytest

from planetgen import tuning
from planetgen.generation import abundance_tier
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.physics import biogenic_oxygen as bo
from planetgen.physics import hydrosphere, surface_pressure
from planetgen.util import draw
from tests.class_support import draw_planet

AREA_M2 = bo.EARTH_AREA_M2


def _world(**fields):
    """A body with Earth's heat, light and warmth and a small salt sea over 15 percent of it."""
    star = SimpleNamespace(luminosity=3.828e26, temperature=5772.0)
    body = SimpleNamespace(hydrosphere="brine seep", ocean_fraction=0.15, surface_temperature=288.0, star=star, distance=1.0,
                           internal_heat_flux_w_m2=bo.EARTH_HEAT_W / AREA_M2, tidal_flux_w_m2=0.0, mantle_delta_iw=bo.EARTH_MANTLE_DELTA_IW)
    for name, value in fields.items():
        setattr(body, name, value)
    return body


def test_a_quiet_interior_takes_a_third_of_the_oxygen_and_any_other_body_takes_what_it_did():
    plain = bo.reductant_flux(_world(), AREA_M2)
    assert plain == pytest.approx(bo.EARTH_REDUCTANT_MOL_YR, rel=1e-6)
    assert bo.reductant_flux(_world(_quiet_interior=True), AREA_M2) == pytest.approx(plain * tuning.ABUNDANCE_QUIET_SINK_FACTOR)
    assert bo.reductant_flux(_world(_quiet_interior=False), AREA_M2) == plain


def test_a_small_sea_makes_its_areas_share_of_an_oceans_oxygen_and_a_seep_makes_a_refuges_small_share():
    sea = bo.productivity(_world(hydrosphere="surface ocean"), AREA_M2)
    assert sea == pytest.approx(bo.EARTH_BURIAL_MOL_YR * 0.15 / bo.SOURCE_REFERENCE_COVER, rel=0.02)
    earth_sea = bo.productivity(_world(hydrosphere="surface ocean", ocean_fraction=surface_pressure.OCEAN_COVER_EARTH), AREA_M2)
    assert earth_sea == pytest.approx(bo.EARTH_BURIAL_MOL_YR, rel=0.02)         # GEN.1837: the generator's Earth, with its own 0.55 of sea, makes Earth's source (a 0.71 reference left it a fifth short)
    seep = bo.productivity(_world(), AREA_M2)
    assert seep == bo.refuge_productivity(_world(), AREA_M2) < 0.1 * sea                # a brine seep is no living sea: only the refuge of GEN.1723, a few percent of the sea's
    assert bo.productivity(_world(_water_dialed="dry"), AREA_M2) == seep                # and the oxygen law does not look at the dial: the hydrosphere law makes the desert sea


def test_with_a_quiet_interior_a_warm_desert_sea_out_produces_its_volcanoes_and_without_one_it_does_not():
    warm = dict(hydrosphere="surface ocean", surface_temperature=300.0)
    source = bo.productivity(_world(**warm), AREA_M2)
    assert bo.reductant_flux(_world(**warm), AREA_M2) < source * 1.2          # Earth's volcanoes take about what a 15 percent sea makes
    assert bo.reductant_flux(_world(_quiet_interior=True, **warm), AREA_M2) < 0.5 * source


def _planets(tier, count, quiet=None, monkeypatch=None):
    if quiet is not None:
        monkeypatch.setitem(tuning.ABUNDANCE_TIER_DIALS, tier, (0.0, 0.0, 0.0, 0.0, quiet, 0.0))
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    config.ABUNDANCE_TIER = tier
    with draw.bound("tier-desert-star"):
        star = Star(config)
    return [draw_planet(star, star.habitable_zone[0] * 1.1, seed, config, moon_count=0) for seed in range(count)]


def test_at_a_quiet_share_of_one_every_candidate_is_drawn_quiet_and_at_zero_and_on_science_no_planet_holds_the_field(monkeypatch):
    quiet = _planets("abundant", 6, quiet=1.0, monkeypatch=monkeypatch)
    candidates = [planet for planet in quiet if abundance_tier.candidate(planet)]
    assert candidates and all(abundance_tier.quiet_interior(planet) for planet in quiet if "_quiet_interior" in vars(planet))
    assert sum("_quiet_interior" in vars(planet) for planet in quiet) >= len(candidates) - 1      # drawn on the zone held at the draw; a planet moved across the edge later keeps the decision
    for tier, share in (("abundant", 0.0), (None, None), ("science", None)):
        planets = _planets(tier, 3, quiet=share, monkeypatch=monkeypatch) if share is not None else _planets(tier, 3)
        assert not any("_quiet_interior" in vars(planet) for planet in planets)          # the golden digests list every field: science holds none
