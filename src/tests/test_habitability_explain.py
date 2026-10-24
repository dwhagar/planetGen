"""
The composition behind the Radiation rating (GEN.445, `planetgen.physics.habitability_explain`): the terms name the air, the field, the
ground's potassium and the ozone, and changing any of them moves the terms and the rating as the model says. Earth is a test target,
never an input of the code.
"""

import types

import pytest

from planetgen.describe import build
from planetgen.physics import (
    body_radiation,
    constants,
    habitability,
    habitability_explain,
    habitability_world,
    radiation,
    radiogenic,
)
from tests.air_support import fold_air

AGE_GY = 4.568
FELSIC = 0.3


def _sun():
    return types.SimpleNamespace(temperature=5772.0, luminosity=constants.SOLAR_LUMINOSITY, uv_c_share=1.0)


def _earth(**fields):
    """An Earth-like record whose dose parts are the model's own (the cosmic-ray curve, the ground of its stored K, Th and U, the
    stored ultraviolet of its air), so a change of the record is carried through the same functions the generator uses."""
    k, th, u = radiogenic.BULK_SILICATE_EARTH
    base = {"planet_class": "M", "body_type": "r", "gravity": 1.0, "atmospheric_pressure": 101325.0, "distance": 1.0, "is_moon": False,
            "star": _sun(), "p_n2_kpa": 79.1, "p_o2_kpa": 21.2, "p_ar_kpa": 0.95, "ozone_column_du": 300.0,
            "surface_temperature": 288.0, "radius": 6371.0, "mass": 5.972e24, "hydrosphere": "ocean",
            "magnetic_moment_a_m2": 7.768e22, "c_k0_ppm": k, "c_th0_ppm": th * radiogenic.formation_mix("Th")[1],
            "c_u0_ppm": u * radiogenic.formation_mix("U")[1], "al26_ratio": 0.0, "urey_ratio": 0.5, "dose_helio_mult": 1.0}
    base.update(fields)
    body = types.SimpleNamespace(**fold_air(base))
    _redo_dose(body)
    return body


def _redo_dose(body):
    """The dose parts and the ultraviolet from the record, as `radiation.generate_dose` draws them with the draws fixed."""
    column = radiation.column_g_cm2(body.atmospheric_pressure, body.gravity)
    cutoff = radiation.cutoff_rigidity_gv(body.magnetic_moment_a_m2, body.radius)
    content = radiogenic.content_of(body)
    scale_height_km = 8.5 if body.atmospheric_pressure else None
    body.dose_gcr_msv_yr = radiation.gcr_dose_msv_yr(column, cutoff, 1.0)
    body.dose_sep_msv_yr = radiation.sep_dose_msv_yr(radiation.free_sep_msv_yr(None, 1.0), column, cutoff)
    body.dose_ground_msv_yr = radiation.ground_dose_msv_yr(content, AGE_GY, FELSIC, scale_height_km=scale_height_km)
    body.surface_dose_msv_yr = body.dose_gcr_msv_yr + body.dose_sep_msv_yr + body.dose_ground_msv_yr
    body.uv_surface_index = body_radiation.dna_dose_index(body)


def _tier_index(terms):
    return ("Blue", "Green", "Yellow", "Red").index(terms.tier)


def test_the_terms_name_the_air_the_field_the_ground_and_the_ozone():
    terms = habitability_explain.radiation_terms(_earth())
    assert terms.column_g_cm2 == pytest.approx(1033.0, rel=0.02)
    assert [name for name, _share in terms.air] == ["nitrogen", "oxygen"]
    assert terms.air[0][1] == pytest.approx(0.78, abs=0.01)
    assert terms.cutoff_gv == pytest.approx(14.34, rel=0.01)
    assert terms.potassium_ppm == pytest.approx(radiogenic.BULK_SILICATE_EARTH[0])
    assert terms.ozone_du == 300.0
    assert terms.total_msv_yr == pytest.approx(terms.gcr_msv_yr + terms.sep_msv_yr + terms.ground_msv_yr)
    assert terms.tier == "Blue"
    assert terms.uv_index > 0.0


def test_a_body_with_no_stored_dose_or_a_giant_has_no_terms():
    with pytest.raises(ValueError):
        habitability_explain.radiation_terms(_earth(surface_dose_msv_yr=None, body_type="g"))
    bare = _earth()
    bare.surface_dose_msv_yr = None
    with pytest.raises(ValueError):
        habitability_explain.radiation_terms(bare)


def test_more_potassium_in_the_ground_moves_the_ground_term_and_the_dose_and_never_the_cosmic_rays():
    base = habitability_explain.radiation_terms(_earth())
    k = radiogenic.BULK_SILICATE_EARTH[0]
    rich = habitability_explain.radiation_terms(_earth(c_k0_ppm=10.0 * k))
    assert rich.potassium_ppm == 10.0 * k
    assert rich.ground_msv_yr > 2.0 * base.ground_msv_yr
    assert rich.total_msv_yr > base.total_msv_yr
    assert rich.gcr_msv_yr == base.gcr_msv_yr and rich.sep_msv_yr == base.sep_msv_yr
    # the gamma law is linear in the element: ten times the potassium adds ten times its own part
    def only_k(ppm):
        return radiation.gamma_dose_msv_yr(radiogenic.surface_content(radiogenic.Content(ppm, 0.0, 0.0, 0.0), FELSIC), AGE_GY)

    assert rich.ground_msv_yr - base.ground_msv_yr == pytest.approx(only_k(10.0 * k) - only_k(k), rel=1e-6)
    # enough potassium turns the dose rating down
    hot = habitability_explain.radiation_terms(_earth(c_k0_ppm=1.0e6))
    assert _tier_index(hot) > _tier_index(base)
    assert hot.limit == "dose"


def test_taking_the_oxygen_and_its_ozone_away_moves_the_ultraviolet_and_the_rating():
    base = habitability_explain.radiation_terms(_earth())
    no_o2 = habitability_explain.radiation_terms(_earth(p_o2_kpa=0.0, p_n2_kpa=100.3, ozone_column_du=0.0))
    assert no_o2.ozone_du == 0.0
    assert no_o2.uv_index > 100.0 * base.uv_index
    assert no_o2.limit == "ultraviolet"
    assert _tier_index(no_o2) > _tier_index(base)
    assert [name for name, _share in no_o2.air] == ["nitrogen"]
    assert no_o2.gcr_msv_yr == pytest.approx(base.gcr_msv_yr, rel=0.01)      # the same column of air stops the same cosmic rays


def test_thinning_the_air_raises_the_cosmic_rays_and_the_stellar_particles():
    base = habitability_explain.radiation_terms(_earth())
    thin = habitability_explain.radiation_terms(_earth(atmospheric_pressure=1013.25, p_n2_kpa=0.79, p_o2_kpa=0.21, p_ar_kpa=0.0095))
    assert thin.column_g_cm2 == pytest.approx(base.column_g_cm2 / 100.0, rel=1e-6)
    assert thin.gcr_msv_yr > 10.0 * base.gcr_msv_yr
    airless = _earth(atmospheric_pressure=0.0, p_n2_kpa=0.0, p_o2_kpa=0.0, p_ar_kpa=0.0)
    terms = habitability_explain.radiation_terms(airless)
    assert terms.column_g_cm2 == 0.0 and terms.air == () and terms.gcr_msv_yr > thin.gcr_msv_yr


def test_a_switched_off_field_raises_the_cosmic_rays():
    base = habitability_explain.radiation_terms(_earth())
    bare = habitability_explain.radiation_terms(_earth(magnetic_moment_a_m2=0.0))
    assert bare.cutoff_gv == 0.0 and base.cutoff_gv > 10.0
    assert bare.gcr_msv_yr > base.gcr_msv_yr


def test_the_rating_follows_the_terms_through_the_same_scores():
    for body in (_earth(), _earth(c_k0_ppm=1.0e6), _earth(p_o2_kpa=0.0, p_n2_kpa=100.3, ozone_column_du=0.0)):
        terms = habitability_explain.radiation_terms(body)
        world = habitability_world.world_from(habitability_world.values_from_body(body))
        assert terms.tier == habitability.tier(habitability.radiation_score(world))
        uv = habitability.uv_score(world)
        assert (terms.limit == "ultraviolet") == (uv is not None and uv < habitability.band_score(world.dose_sv_yr, habitability.HUMAN.dose, log=True))


def test_the_describe_block_writes_every_term_with_its_composition():
    planet = _earth()
    block = build._radiation_block(planet)
    assert block["caption"] == "Radiation"
    keys = [row["key"] for row in block["rows"]]
    assert keys == ["radiation-dose", "radiation-limit", "radiation-cosmic", "radiation-stellar", "radiation-ground", "radiation-uv"]
    text = {row["key"]: "".join(item if isinstance(item, str) else "" for item in row["value"]) for row in block["rows"]}
    assert "nitrogen 78.1 percent" in text["radiation-cosmic"] and "14.3 GV" in text["radiation-cosmic"]
    assert "potassium" in text["radiation-ground"] and "uranium" in text["radiation-ground"]
    assert "300 Dobson units of ozone" in text["radiation-uv"]


def test_the_describe_block_says_so_for_no_air_no_field_and_no_ozone():
    planet = _earth(atmospheric_pressure=0.0, p_n2_kpa=0.0, p_o2_kpa=0.0, p_ar_kpa=0.0, magnetic_moment_a_m2=0.0, ozone_column_du=0.0)
    block = build._radiation_block(planet)
    text = " ".join("".join(item if isinstance(item, str) else "" for item in row["value"]) for row in block["rows"])
    assert "no air to stop them" in text and "no magnetic field" in text


def test_a_gas_giant_or_an_unscored_body_gets_no_block():
    assert build._radiation_block(_earth(body_type="g")) is None
    unscored = _earth()
    unscored.surface_dose_msv_yr = None
    assert build._radiation_block(unscored) is None
