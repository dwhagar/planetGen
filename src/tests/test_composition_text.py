# tests/test_composition_text.py

"""
GEN.430, GEN.642: a body's composition phrase is written from its stored bulk fractions by fixed thresholds, and a giant's page lists its
layers.
"""

from types import SimpleNamespace

import pytest

from planetgen.describe import build, markdown
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import bulk_composition as bc, composition_text, giant_interior
from planetgen.util import draw


def _body(**fractions):
    shares = dict(metal=0.0, rock=0.0, water=0.0, ice=0.0, envelope=0.0)
    shares.update(fractions)
    extra = {key: fractions.pop(key) for key in ("carbon_rock", "porosity") if key in fractions}
    shares = {key: value for key, value in shares.items() if key not in extra}
    return SimpleNamespace(**bc.BulkFractions(**shares, **extra).to_attributes())


def test_the_largest_share_comes_first_and_the_rest_follow_by_level():
    body = _body(metal=0.30, rock=0.62, water=0.06, ice=0.02)
    assert composition_text.text(body) == (
        "mostly silicate rock, with a large share of iron-rich metal, a moderate share of water and a small share of volatile ices")


def test_materials_of_the_same_level_are_named_together():
    body = _body(metal=0.10, rock=0.40, water=0.40, ice=0.10)
    assert composition_text.text(body) == "a large share of silicate rock and water, with a moderate share of iron-rich metal and volatile ices"


def test_a_material_under_one_percent_is_not_named():
    body = _body(metal=0.2, rock=0.795, water=0.005)
    assert "water" not in composition_text.text(body)


def test_a_carbon_rich_body_names_silicon_carbide_and_graphite():
    assert composition_text.text(_body(metal=0.3, rock=0.7, carbon_rock=True)).startswith("mostly silicon carbide and graphite")


def test_a_porous_small_body_is_loosely_packed():
    assert composition_text.text(_body(rock=0.5, water=0.5, porosity=0.2)).endswith(", loosely packed")
    assert not composition_text.text(_body(rock=0.5, water=0.5, porosity=0.02)).endswith("packed")


def test_a_body_with_no_stored_composition_has_no_text():
    assert composition_text.text(SimpleNamespace()) is None


def test_every_generated_body_states_what_its_fractions_say():
    seen = set()
    for seed in range(6):
        with draw.bound(seed):
            system = StarSystem(system_config=SystemConfig())
        for planet in system.planets:
            if not hasattr(planet, "moons"):
                continue
            for body in [planet, *planet.moons]:
                bulk = bc.bulk_fractions(body)
                assert body.composition == composition_text.text(body)
                for fraction, share in bulk.shares().items():
                    name = composition_text.rock_name(bulk.carbon_rock) if fraction == "rock" else composition_text.NAMES[fraction]
                    assert (name in body.composition) == (share >= composition_text.MIN_SHARE), (body.planet_class, fraction, share)
                seen.add(body.body_type)
    assert seen == {"t", "g"}


def test_a_giants_page_lists_its_layers_and_a_rocky_worlds_does_not():
    for seed in range(60):
        with draw.bound(seed):
            system = StarSystem(system_config=SystemConfig())
        giants = [p for p in system.planets if getattr(p, "body_type", "") == "g"]
        if giants:
            break
    else:
        pytest.fail("no giant in 60 seeds")
    page = markdown.render(build.describe(system))
    for sentence in giant_interior.describe_interior(giant_interior.interior_of(giants[0])):
        assert sentence in page
    assert "bar down to" in page or "GPa down to" in page
