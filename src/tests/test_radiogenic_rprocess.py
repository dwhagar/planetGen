# tests/test_radiogenic_rprocess.py

"""
GEN.759: the thorium to uranium ratio of a body's rock follows the r-process enrichment of its host's birth gas
(`physics.rprocess_enrichment`, GEN.663) instead of a fixed scatter on the bulk silicate Earth's uranium.
"""

import math
import pathlib
from types import SimpleNamespace

import pytest

from planetgen.galaxy import shape_draw
from planetgen.generation.config import SystemConfig
from planetgen.physics import radiogenic, rprocess_enrichment
from planetgen.util import draw

SOURCE = pathlib.Path(radiogenic.__file__).parent


def test_the_ratio_the_formula_gives_brackets_the_bulk_silicate_earth():
    low, high = rprocess_enrichment.PRODUCTION_RATIO_RANGE
    factor = radiogenic.default_th_u_factor()
    assert 1.0 < factor < 1.6
    bse = radiogenic.BULK_SILICATE_EARTH[1] / radiogenic.BULK_SILICATE_EARTH[2]
    assert radiogenic.th_u_ratio_now(low, factor) < bse < radiogenic.th_u_ratio_now(high, factor) * 1.01   # a test point, not an input


def test_the_ratio_is_the_published_formula_at_the_solar_systems_age():
    thorium, uranium = rprocess_enrichment.decay_constants_per_gy()
    decay = math.exp((uranium - thorium) * radiogenic.REFERENCE_AGE_GY)
    # 0.97483 is the mass ratio of the two isotopes; the element's uranium is 99.27 percent 238U today
    expected = 0.97483 * 1.8 * 1.2 * decay / 0.9927
    assert radiogenic.th_u_ratio_now(1.8, 1.2) == pytest.approx(expected, rel=2e-3)


def test_a_body_draws_the_production_ratio_where_it_drew_the_uranium_scatter_and_keeps_the_stream_in_step():
    a, b = draw.Stream("rp-a"), draw.Stream("rp-a")
    contents = [radiogenic.draw_content(0.0, rng=a, th_u_factor=1.2) for _ in range(200)]
    for _ in range(200):
        radiogenic.draw_content(0.3, rng=b, th_u_factor=1.4)
    assert a.random() == b.random()
    ratios = [c.th0_ppm / c.u0_ppm * radiogenic.formation_mix("U")[1] / radiogenic.formation_mix("Th")[1] for c in contents]
    low, high = rprocess_enrichment.PRODUCTION_RATIO_RANGE
    assert min(ratios) == pytest.approx(radiogenic.th_u_ratio_now(low, 1.2), rel=0.02)
    assert max(ratios) == pytest.approx(radiogenic.th_u_ratio_now(high, 1.2), rel=0.02)


def test_gas_enriched_long_before_the_star_formed_leaves_the_rock_thorium_rich_and_uranium_poor():
    same = dict(rng=draw.Stream("rp-b"))
    poor = radiogenic.draw_content(0.0, th_u_factor=1.0, **same)
    rich = radiogenic.draw_content(0.0, th_u_factor=1.5, rng=draw.Stream("rp-b"))
    assert rich.th0_ppm == pytest.approx(poor.th0_ppm)             # thorium follows the star's Th/Mg
    assert rich.u0_ppm == pytest.approx(poor.u0_ppm / 1.5)         # the uranium is what the ratio leaves
    assert rich.k0_ppm == pytest.approx(poor.k0_ppm / 1.5)         # potassium is a multiple of the uranium


def test_a_younger_host_in_the_same_galaxy_has_the_larger_factor():
    shape = shape_draw.default_shape(None)
    # a star formed later finds more earlier ejecta already decayed, so its gas is thorium-richer
    factors = [rprocess_enrichment.th_u_factor_of_shape(shape, age) for age in (9.0, 5.0, 1.0)]
    assert factors[0] < factors[1] < factors[2] and factors[0] >= 1.0


def test_a_body_reads_the_factor_of_its_hosts_age_in_its_galaxy():
    shape = shape_draw.default_shape(None)
    config = SystemConfig()
    config.GALAXY_SHAPE = shape
    body = SimpleNamespace(star=SimpleNamespace(age=7.0), system_config=config)
    assert radiogenic.star_th_u_factor(body) == rprocess_enrichment.th_u_factor_of_shape(shape, 7.0)
    config.GALAXY_SHAPE = None                                       # no stored galaxy: the default shape, as any stand-alone work takes
    assert radiogenic.star_th_u_factor(body) == rprocess_enrichment.th_u_factor_of_shape(shape, 7.0)
    assert radiogenic.star_th_u_factor(SimpleNamespace(star=SimpleNamespace(age=None), system_config=config)) is None


def test_the_typical_star_lies_between_the_newest_and_the_oldest():
    shape = shape_draw.default_shape(None)
    typical = rprocess_enrichment.typical_th_u_factor(shape)
    assert rprocess_enrichment.th_u_factor_of_shape(shape, 12.0) < typical < rprocess_enrichment.th_u_factor_of_shape(shape, 0.5)


def test_the_enrichment_module_has_a_caller_and_the_old_scatter_is_gone():
    text = (SOURCE / "radiogenic.py").read_text(encoding="utf-8")
    assert "rprocess_enrichment." in text and "th_u_factor_of_shape" in text            # the wiring guard (GEN.754)
    tuning_text = (SOURCE.parent / "tuning.py").read_text(encoding="utf-8")
    assert "RADIOGENIC_URANIUM_EXTRA_SIGMA_DEX" not in tuning_text and "uranium_extra" not in text
