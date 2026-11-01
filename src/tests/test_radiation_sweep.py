"""
A radiation sweep over every class (TEST.212 part 1): band transmission, greenhouse, shielding and dose agree with the composition.

Generates systems in memory (the census harness's `sample_bodies`, half of them asked for a habitable world so the rare temperate classes appear) and
checks every planet and moon the generator stored against the rules the radiation code states. Each rule is a function that returns the problems it
finds in one body (`*_problems`), so the sweep asserts there are none and a second set of tests breaks each rule in turn on a copy of a good body and
asserts the rule names the break.

The rules:

1. **Transmission.** Every band of the light at the ground (`body_radiation.surface_spectrum`) is between 0 and 1 of the light at the top of the air; a
   body with air and a lit star has a spectrum and one without has none; a thicker column of the same air passes no more light than the thinner one
   (the transmission follows the mixture).
2. **Greenhouse.** The stored Bond albedo is between 0 and 1; the greenhouse factor is zero for a body with no air and positive for one with air; the
   stored surface temperature is `T_e (1 + G)^(1/4)` of the stored albedo, greenhouse factor and starlight (`climate_draws.surface_temperature_k`), or,
   for a surface held at `surface_air.MAX_SURFACE_TEMPERATURE_K`, no hotter than that; and over the rocky planets with air the greenhouse rises with the
   infrared-active column (`gases.Gas.is_infrared_active`).
3. **Dose.** The yearly dose is the sum of its four terms; the cosmic-ray term lies between the dose the column and the field give at the two ends of
   the solar-modulation range; the stellar-particle term is the one the column, the field and the star's flares give; a giant has no dose; the belt dose
   is a moon's only; the ground dose of a rocky body is positive.
4. **Ultraviolet.** The stored index of a planet is the Setlow dose of the spectrum below its air as a multiple of Earth's (`dna_dose_index`); the
   ozone column is positive exactly when the air holds oxygen under a lit star; a body with air takes no more ultraviolet than the bare star's light.

A moon whose planet was moved after the moon's dose was made keeps the dose of the old distance (reported to TODO), so the rules that read the distance
are applied to planets and the others to every body.

Not here (they wait for GEN.444's later slices and GEN.446): the Earth, Moon and Mars doses against the measured values, and the descriptions that
quote these numbers.
"""

import collections
import copy
import math

import pytest
from scipy import stats

from planetgen import tuning
from planetgen.physics import (air_column, atmosphere, body_radiation, climate_draws, constants, gases, radiation as rad, surface_air)
from tests.class_sweep_cli import sample_bodies

SYSTEMS = 14
"""int: Systems generated for the sweep (about 90 s; 160 or so bodies of 15 classes)."""
SEED = 11
HABITABLE_SHARE = 0.5
RELATIVE_TOLERANCE = 1e-6
"""float: What counts as the same stored number."""
BANDS = ("uv_c", "uv_b", "uv_a", "visible", "near_infrared")
MIN_CLASSES = 8
"""int: The classes the sweep must reach, or it has not swept every kind of body."""


@pytest.fixture(scope="module")
def swept():
    """`(body, is_moon)` for every planet and moon of the sweep."""
    return sample_bodies(SYSTEMS, SEED, HABITABLE_SHARE)


def _starlight_w_m2(body):
    """The starlight at the body, W/m2, or `None` for a body with no star light or no stored distance."""
    distance_au = getattr(body, "_star_distance_au", None)
    luminosity = getattr(getattr(body, "star", None), "luminosity", None)
    if not distance_au or not luminosity:
        return None
    return climate_draws.starlight_w_m2(luminosity, distance_au * constants.AU_TO_KM)


def _same(a, b):
    return math.isclose(a, b, rel_tol=RELATIVE_TOLERANCE, abs_tol=1e-12)


def _has_air(body):
    return (getattr(body, "atmospheric_pressure", None) or 0.0) > 0.0


def transmission_problems(body):
    """Rule 1 for one body: the band transmissions in 0..1, the spectrum present exactly when there is air and starlight, and a thicker column of
    the same air passing no more of the light."""
    spectrum = body_radiation.surface_spectrum(body)
    lit = _starlight_w_m2(body) is not None
    if spectrum is None:
        return ["has air and a lit star but no spectrum"] if (air_column.column_of(body) is not None and lit) else []
    problems = []
    for band, value in spectrum.band_fractions().items():
        if not 0.0 <= value <= 1.0 + 1e-6:
            problems.append(f"band {band} transmission {value:.4g} is outside 0..1")
    thicker = copy.copy(body)
    thicker.atmospheric_pressure = body.atmospheric_pressure * 2.0
    setattr(thicker, atmosphere.GAS_FIELD, {formula: pa * 2.0 for formula, pa in atmosphere.gas_partials_pa(body).items()})
    other = body_radiation.surface_spectrum(thicker)
    if other is not None:
        for band in ("uv_c", "uv_b", "uv_a", "visible"):
            if other.band_transmission(band) > spectrum.band_transmission(band) * (1.0 + 1e-3) + 1e-9:
                problems.append(f"a column twice as thick passes more {band} light")
    return problems


def climate_problems(body):
    """Rule 2 for one body, except the infrared ranking, which needs the whole sweep (`test_the_greenhouse_rises_with_the_infrared_active_column`)."""
    problems = []
    albedo, factor, temperature = body.bond_albedo, body.greenhouse_factor, body.surface_temperature
    if albedo is None or not 0.0 < albedo < 1.0:
        problems.append(f"Bond albedo {albedo} is outside 0..1")
        return problems
    if factor is None or not math.isfinite(factor) or factor < 0.0:
        problems.append(f"greenhouse factor {factor} is not a finite number of zero or more")
        return problems
    if not _has_air(body) and factor != 0.0:
        problems.append("an airless body has a greenhouse factor")
    flux = _starlight_w_m2(body)
    if flux is not None:
        expected = climate_draws.surface_temperature_k(albedo, factor, flux)
        cap = surface_air.MAX_SURFACE_TEMPERATURE_K
        if _same(temperature, cap):
            if expected < cap * (1.0 - RELATIVE_TOLERANCE):
                problems.append(f"surface held at the {cap} K cap although albedo and greenhouse give {expected:.1f} K")
        elif not math.isclose(temperature, expected, rel_tol=1e-4):
            problems.append(f"stored {temperature:.2f} K but albedo and greenhouse give {expected:.2f} K")
    return problems


def dose_problems(body, is_moon):
    """Rule 3 for one body."""
    fields = (body.surface_dose_msv_yr, body.dose_gcr_msv_yr, body.dose_sep_msv_yr, body.dose_ground_msv_yr, body.dose_belt_msv_yr)
    if body.body_type == "g":
        return [] if all(value is None for value in fields) else ["a giant has a surface dose"]
    problems = []
    total, gcr, sep, ground, belt = fields
    if not _same(total, gcr + sep + ground + belt):
        problems.append(f"dose {total:.6g} is not the sum of its terms {gcr + sep + ground + belt:.6g}")
    if min(gcr, sep, ground, belt) < 0.0 or gcr <= 0.0 or ground <= 0.0:
        problems.append("a dose term is negative, or the cosmic-ray or ground term is not positive")
    if belt > 0.0 and not is_moon:
        problems.append("a planet has a belt dose")
    if is_moon:
        return problems
    column = rad.column_g_cm2(getattr(body, "atmospheric_pressure", None) or 0.0, getattr(body, "gravity", None))
    cutoff = rad.cutoff_rigidity_gv(getattr(body, "magnetic_moment_a_m2", None) or 0.0, body.radius)
    low, high = (rad.gcr_dose_msv_yr(column, cutoff, share) for share in tuning.GCR_SOLAR_MODULATION_RANGE)
    if not low * (1.0 - 1e-9) <= gcr <= high * (1.0 + 1e-9):
        problems.append(f"cosmic-ray dose {gcr:.4g} is outside {low:.4g} to {high:.4g}, what the column and the field give")
    star = body.star
    free = rad.free_sep_msv_yr(getattr(star, "flare_n33_per_yr", None), getattr(body, "_star_distance_au", None), getattr(star, "flare_alpha", None))
    if not math.isclose(sep, rad.sep_dose_msv_yr(free, column, cutoff), rel_tol=1e-6, abs_tol=1e-12):
        problems.append(f"stellar-particle dose {sep:.4g} is not what the column, the field and the flares give")
    return problems


def ultraviolet_problems(body, is_moon):
    """Rule 4 for one body (not a giant)."""
    problems = []
    has_light = bool(getattr(body.star, "temperature", None))
    holds_oxygen = atmosphere.partial_kpa(body, "o2") > 0.0
    if (body.ozone_column_du > 0.0) != (holds_oxygen and has_light):
        problems.append("the ozone column is not positive exactly when the air holds oxygen under a lit star")
    index = body_radiation.dna_dose_index(body)
    if not is_moon and index is not None and not _same(body.uv_surface_index, index):
        problems.append(f"stored UV index {body.uv_surface_index:.6g} is not the Setlow dose of the spectrum {index:.6g}")
    bare = body_radiation.bare_dna_dose_w_m2(body)
    if _has_air(body) and bare is not None and not is_moon:
        spectrum = body_radiation.surface_spectrum(body)
        if spectrum is not None and spectrum.dna_dose_w_m2() > bare * (1.0 + 1e-6):
            problems.append("a body with air takes more ultraviolet than the bare star's light")
    return problems


def test_the_sweep_reaches_every_kind_of_body(swept):
    classes = collections.Counter(body.planet_class for body, _ in swept)
    kinds = {body.body_type for body, _ in swept}
    assert len([letter for letter in classes if letter]) >= MIN_CLASSES, classes
    assert {"t", "g"} <= kinds
    assert any(is_moon for _, is_moon in swept) and any(not is_moon for _, is_moon in swept)
    assert any(_has_air(body) for body, _ in swept) and any(not _has_air(body) for body, _ in swept)


def test_every_band_transmission_follows_the_mixture(swept):
    problems = {(body.planet_class, text) for body, _ in swept for text in transmission_problems(body)}
    assert not problems, sorted(problems, key=str)


def test_every_surface_temperature_is_the_albedo_and_greenhouse_the_starlight_gives(swept):
    problems = {(body.planet_class, text) for body, _ in swept for text in climate_problems(body)}
    assert not problems, sorted(problems, key=str)


def test_the_greenhouse_rises_with_the_infrared_active_column(swept):
    columns, factors = [], []
    for body, _ in swept:
        column = air_column.column_of(body)
        if body.body_type != "t" or column is None:
            continue
        active = sum(column.number_column_cm2(formula) for formula in column.partials_pa if gases.gas(formula).is_infrared_active)
        columns.append(math.log10(1.0 + active))
        factors.append(body.greenhouse_factor)
    assert len(columns) >= 20
    rho = stats.spearmanr(columns, factors).statistic
    assert rho > 0.3, f"the greenhouse factor follows the infrared-active column with rank correlation {rho:.2f}"


def test_every_dose_is_its_terms_and_follows_the_column_the_field_and_the_flares(swept):
    problems = {(body.planet_class, text) for body, is_moon in swept for text in dose_problems(body, is_moon)}
    assert not problems, sorted(problems, key=str)


def test_every_ultraviolet_index_is_the_dose_of_the_light_below_the_air(swept):
    problems = {(body.planet_class, text) for body, is_moon in swept if body.body_type != "g" for text in ultraviolet_problems(body, is_moon)}
    assert not problems, sorted(problems, key=str)


def _first(swept, test):
    for body, is_moon in swept:
        if test(body, is_moon):
            return copy.copy(body), is_moon
    pytest.fail("the sweep holds no body to break")


def test_a_transmission_outside_zero_to_one_is_named(swept, monkeypatch):
    body, _ = _first(swept, lambda b, m: air_column.column_of(b) is not None and _starlight_w_m2(b) is not None)
    spectrum = body_radiation.surface_spectrum(body)
    monkeypatch.setattr(type(spectrum), "band_fractions", lambda self: {band: 1.5 for band in BANDS})
    assert any("outside 0..1" in text for text in transmission_problems(body))


def test_a_missing_spectrum_under_a_lit_star_is_named(swept, monkeypatch):
    body, _ = _first(swept, lambda b, m: air_column.column_of(b) is not None and _starlight_w_m2(b) is not None)
    monkeypatch.setattr(body_radiation, "surface_spectrum", lambda *args, **kwargs: None)
    assert any("no spectrum" in text for text in transmission_problems(body))


def test_a_thicker_column_that_passes_more_light_is_named(swept, monkeypatch):
    body, _ = _first(swept, lambda b, m: air_column.column_of(b) is not None and _starlight_w_m2(b) is not None and b.atmospheric_pressure > 1e3)
    real = body_radiation.surface_spectrum(body)

    def fewer_photons_for_thin_air(candidate, *args, **kwargs):
        return real if candidate.atmospheric_pressure == body.atmospheric_pressure else _Brighter(real)

    monkeypatch.setattr(body_radiation, "surface_spectrum", fewer_photons_for_thin_air)
    assert any("passes more" in text for text in transmission_problems(body))


class _Brighter:
    """A spectrum that passes ten times the light of `real` in every band."""

    def __init__(self, real):
        self._real = real

    def band_fractions(self):
        return {band: value * 10.0 for band, value in self._real.band_fractions().items()}

    def band_transmission(self, band):
        return 10.0 * self._real.band_transmission(band)


def test_a_greenhouse_that_disagrees_with_the_temperature_is_named(swept):
    body, _ = _first(swept, lambda b, m: _has_air(b) and _starlight_w_m2(b) is not None and b.surface_temperature < 0.9 * surface_air.MAX_SURFACE_TEMPERATURE_K)
    body.greenhouse_factor += 1.0
    assert any("albedo and greenhouse give" in text for text in climate_problems(body))


def test_a_surface_held_at_the_cap_with_a_cooler_greenhouse_is_named(swept):
    body, _ = _first(swept, lambda b, m: _has_air(b) and _starlight_w_m2(b) is not None and b.surface_temperature < 0.5 * surface_air.MAX_SURFACE_TEMPERATURE_K)
    body.surface_temperature = surface_air.MAX_SURFACE_TEMPERATURE_K
    assert any("cap" in text for text in climate_problems(body))


def test_an_albedo_outside_zero_to_one_and_a_greenhouse_over_no_air_are_named(swept):
    body, _ = _first(swept, lambda b, m: not _has_air(b))
    body.greenhouse_factor = 3.0
    assert any("airless" in text for text in climate_problems(body))
    body.bond_albedo = 1.2
    assert any("Bond albedo" in text for text in climate_problems(body))


def test_a_dose_that_is_not_its_terms_is_named(swept):
    body, is_moon = _first(swept, lambda b, m: b.body_type != "g")
    body.surface_dose_msv_yr *= 2.0
    assert any("not the sum" in text for text in dose_problems(body, is_moon))


def test_a_cosmic_ray_dose_the_column_does_not_give_is_named(swept):
    body, is_moon = _first(swept, lambda b, m: b.body_type != "g" and not m)
    body.dose_gcr_msv_yr *= 100.0
    body.surface_dose_msv_yr = body.dose_gcr_msv_yr + body.dose_sep_msv_yr + body.dose_ground_msv_yr + body.dose_belt_msv_yr
    assert any("cosmic-ray dose" in text for text in dose_problems(body, is_moon))


def test_a_stellar_particle_dose_the_flares_do_not_give_is_named(swept):
    body, is_moon = _first(swept, lambda b, m: b.body_type != "g" and not m and b.dose_sep_msv_yr > 0.0)
    body.dose_sep_msv_yr *= 3.0
    body.surface_dose_msv_yr = body.dose_gcr_msv_yr + body.dose_sep_msv_yr + body.dose_ground_msv_yr + body.dose_belt_msv_yr
    assert any("stellar-particle dose" in text for text in dose_problems(body, is_moon))


def test_a_giant_with_a_dose_and_a_planet_with_a_belt_are_named(swept):
    giant, is_moon = _first(swept, lambda b, m: b.body_type == "g")
    giant.surface_dose_msv_yr = 1.0
    assert any("giant" in text for text in dose_problems(giant, is_moon))
    planet, is_moon = _first(swept, lambda b, m: b.body_type != "g" and not m)
    planet.dose_belt_msv_yr = 5.0
    planet.surface_dose_msv_yr += 5.0
    assert any("belt" in text for text in dose_problems(planet, is_moon))


def test_an_ultraviolet_index_that_is_not_the_spectrums_is_named(swept):
    body, is_moon = _first(swept, lambda b, m: b.body_type != "g" and not m and body_radiation.dna_dose_index(b) is not None)
    body.uv_surface_index = body.uv_surface_index * 2.0 + 1.0
    assert any("UV index" in text for text in ultraviolet_problems(body, is_moon))


def test_an_ozone_column_with_no_oxygen_is_named(swept):
    body, is_moon = _first(swept, lambda b, m: b.body_type != "g" and atmosphere.partial_kpa(b, "o2") == 0.0 and getattr(b.star, "temperature", None))
    body.ozone_column_du = 100.0
    assert any("ozone column" in text for text in ultraviolet_problems(body, is_moon))
