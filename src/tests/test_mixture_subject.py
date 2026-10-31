"""
The description of a body the atmosphere mixture draw reads (GEN.756, `planetgen.physics.mixture_subject`). The Earth is a test
target, never an input of the code.
"""

import math
from types import SimpleNamespace

import pytest

from planetgen import tuning
from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.star import Star
from planetgen.physics import activity, constants, mixture_draw, mixture_subject as ms, thermosphere
from planetgen.physics.mixture_base import MixtureError
from planetgen.util import draw

EARTH_MASS_KG, EARTH_RADIUS_KM = constants.EARTH_MASS_TO_KG, constants.EARTH_RADIUS_KM


def _sun(**over):
    values = dict(type="G2V", yerkes_class="V", luminosity=constants.SOLAR_LUMINOSITY, mass=constants.SOLAR_MASS_TO_KG, temperature=5772.0,
                  age=4.5, habitable_zone=(0.95, 1.68), snow_line_au=2.7, xuv_fluence_j=activity.REFERENCE_FLUENCE_J)
    values.update(over)
    return SimpleNamespace(**values)


def _earth(**over):
    values = dict(body_type="t", mass=EARTH_MASS_KG, radius=EARTH_RADIUS_KM, distance=1.0, bond_albedo=0.306, is_moon=False,
                  star=_sun(), core_mass_earth=None)
    values.update(over)
    return SimpleNamespace(**values)


def test_the_earth_gives_the_quantities_the_world_gives():
    q = ms.quantities_of_body(_earth(), ts_k=288.0, pressure_bar=1.0, ocean_present=True, long_lived_geology=True)
    assert q["mass_earth"] == pytest.approx(1.0)
    assert q["v_esc_kms"] == pytest.approx(11.19, rel=0.01)
    assert q["insolation"] == pytest.approx(1.0)
    assert q["teq_k"] == pytest.approx(254.0, abs=2.0)
    assert q["xuv_index"] == pytest.approx(1.0)
    assert q["snow_ratio"] == pytest.approx(1.0 / 2.7)
    assert q["zone"] == "e" and q["host_type"] == "G"
    assert q["age_gy"] == 4.5 and q["luminosity_solar"] == pytest.approx(1.0)
    assert q["runaway_ratio"] == pytest.approx(0.90, abs=0.03)
    assert q["tday_k"] == pytest.approx(q["teq_k"] * (8.0 / 3.0) ** 0.25)
    assert q["ocean_present"] is True and q["water_lost"] is False


def test_the_xuv_index_is_the_one_stored_for_the_distance_the_body_sits_at():
    """GEN.806: a body that stored its exposure index at its distance gives that value, and the one implementation answers at another distance."""
    stored = _earth(xuv_exposure_index=7.5, _star_distance_au=1.0)
    assert ms.quantities_of_body(stored, ts_k=288.0, pressure_bar=1.0)["xuv_index"] == 7.5
    assert ms.quantities_of_body(stored, ts_k=288.0, pressure_bar=1.0, star_distance_au=2.0)["xuv_index"] == pytest.approx(0.25)
    unstored = _earth(xuv_exposure_index=None)
    assert ms.quantities_of_body(unstored, ts_k=288.0, pressure_bar=1.0)["xuv_index"] == pytest.approx(1.0)


def test_the_subject_holds_exactly_what_the_class_reads():
    subject = ms.subject_of_body(_earth(), ts_k=288.0, pressure_bar=1.0, ocean_present=True, long_lived_geology=True)
    assert subject.body_class == "terrestrial"
    required = set(mixture_draw.required_inputs("terrestrial"))
    assert set(subject.quantities) == required | {"host_mass_solar"}       # GEN.1199: the escape filter finds the XUV from the star's mass


def test_the_zone_follows_the_habitable_zone():
    assert ms.quantities_of_body(_earth(distance=0.4), ts_k=500.0, pressure_bar=1.0)["zone"] == "h"
    assert ms.quantities_of_body(_earth(distance=1.2), ts_k=288.0, pressure_bar=1.0)["zone"] == "e"
    assert ms.quantities_of_body(_earth(distance=5.0), ts_k=100.0, pressure_bar=1.0)["zone"] == "c"


def test_more_starlight_means_more_insolation_and_a_higher_equilibrium_temperature():
    near = ms.quantities_of_body(_earth(distance=0.5), ts_k=300.0, pressure_bar=1.0)
    far = ms.quantities_of_body(_earth(distance=2.0), ts_k=300.0, pressure_bar=1.0)
    assert near["insolation"] == pytest.approx(16.0 * far["insolation"])
    assert near["teq_k"] > far["teq_k"]
    assert near["runaway_ratio"] > far["runaway_ratio"]


@pytest.mark.parametrize("star,expected", [
    (dict(type="F5V"), "F"), (dict(type="K2V"), "K"), (dict(type="M4V"), "M"), (dict(type="A1V"), "F"), (dict(type="B2V"), "F"),
    (dict(type="L3"), "M"), (dict(type="M5VII", yerkes_class="VII"), "WD"), (dict(type="DA", yerkes_class="D"), "WD"),
])
def test_the_host_type_of_a_star(star, expected):
    assert ms.host_type_of(_sun(**star)) == expected


def test_a_pulsar_is_a_pulsar():
    class NeutronStar:
        type = "Neutron Star (Millisecond Pulsar)"
    assert ms.host_type_of(NeutronStar()) == "PSR"


def test_a_black_hole_takes_the_bare_remnant_row():
    class BlackHole:
        type = "Stellar-Mass Black Hole"
    assert ms.host_type_of(BlackHole()) == "PSR"


def test_a_star_with_no_type_is_refused():
    with pytest.raises(MixtureError):
        ms.host_type_of(_sun(type=None))


def test_the_age_of_a_starless_body_is_not_invented_and_what_needs_the_star_is_refused_by_name():
    rogue = _earth(star=None)
    with pytest.raises(MixtureError) as error:
        ms.subject_of_body(rogue, ts_k=288.0, pressure_bar=1.0)
    assert "age_gy" in str(error.value) or "insolation" in str(error.value)


def test_a_moon_needs_the_distance_of_its_parent_from_the_star():
    moon = _earth(is_moon=True, distance=0.003)
    with pytest.raises(MixtureError) as error:
        ms.subject_of_body(moon, ts_k=288.0, pressure_bar=1.0)
    assert "insolation" in str(error.value)
    subject = ms.subject_of_body(moon, ts_k=288.0, pressure_bar=1.0, star_distance_au=5.2)
    assert subject["insolation"] == pytest.approx(1.0 / 5.2 ** 2)


def test_a_star_with_no_age_takes_the_reference_age():
    q = ms.quantities_of_body(_earth(star=_sun(age=None)), ts_k=288.0, pressure_bar=1.0)
    assert q["age_gy"] == tuning.ROGUE_EARTH_REFERENCE_AGE_GY


def test_a_giant_gets_its_layer_temperatures_and_the_photosphere_pressure():
    jupiter = _earth(body_type="g", mass=317.8 * EARTH_MASS_KG, radius=71492.0, distance=5.2)
    q = ms.quantities_of_body(jupiter, ts_k=165.0, pressure_bar=1.0)
    assert q["mass_jupiter"] == pytest.approx(1.0, rel=0.01)
    assert q["photosphere_bar"] == ms.PHOTOSPHERE_BAR
    assert 100.0 < q["teff_k"] < 200.0
    assert q["t_layer_k"] == pytest.approx(1.33 * q["teff_k"])


@pytest.fixture(scope="module")
def host():
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    return Star(config)


def _real(host, mass_earth, zone):
    inner, outer = host.habitable_zone
    distance = {"h": 0.5 * inner, "e": 0.5 * (inner + outer), "c": 2.0 * outer}[zone]
    return Planet(SystemConfig(), host, host.habitable_zone, distance, mass=mass_earth * constants.EARTH_MASS_TO_KG)


# An Earth-mass rocky planet in the ecosphere, a gas giant and an ice giant in the cold zone.
@pytest.mark.parametrize("mass_earth,zone", [(1.0, "e"), (300.0, "c"), (17.0, "c")])
def test_a_generated_body_gives_a_subject_the_draw_accepts(host, mass_earth, zone):
    planet = _real(host, mass_earth, zone)
    subject = ms.subject_of_body(planet, ts_k=planet.surface_temperature, pressure_bar=max(planet.atmospheric_pressure, 1.0) / 1e5,
                                 long_lived_geology=True)
    required = set(mixture_draw.required_inputs(subject.body_class))
    assert required <= set(subject.quantities) <= required | set(ms.OPTIONAL_INPUTS)       # GEN.804: the star's own activity rides along
    assert subject.quantities["log_lx_lbol"] == host.log_lx_lbol
    assert subject.quantities["host_mass_solar"] == pytest.approx(host.mass / constants.SOLAR_MASS_TO_KG)      # GEN.1199: the escape filter needs the star to find the XUV
    assert thermosphere.xuv_flux_earth(subject) > 0.0
    assert "magnetopause_rp" not in subject.quantities
    with draw.bound(draw.Stream("mixture-subject")):
        result = mixture_draw.draw_mixture(subject)
    assert result.mole_fractions and sum(result.mole_fractions.values()) == pytest.approx(1.0)
    assert all(math.isfinite(float(value)) for value in subject.quantities.values() if not isinstance(value, (str, bool)))
