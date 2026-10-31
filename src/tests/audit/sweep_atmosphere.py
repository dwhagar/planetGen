# tests/audit/sweep_atmosphere.py

"""
The rule set of TEST.209 on the body sweep (`tests.audit.sweep`): the air, the temperature, the ozone, the aerosols and the text of a
planet or moon all agree.

Every rule is written from the laws, not from the generator's own functions where the law is short (the sum of the partial pressures, the
mean molar mass, the ideal gas density and scale height, the Stefan-Boltzmann temperature); the ozone and the aerosol rules read the
physics modules for the bounds only they know (the ozone fit and its spread, the most condensate the air can hold).

Rules (names as `RULES` lists them):

- `pressure-sum`: the partial pressures add up to the stored surface pressure; an airless body has none of either.
- `mean-molar-mass`: the stored mean molar mass is the sum of mole fraction times the gas table's molar mass.
- `air-density`: the stored density is P M / (R T) at the surface.
- `scale-height`: the stored scale height is R T / (M g).
- `temperature`: the stored temperature is ((1 - A) S (1 + G) / (4 sigma))^(1/4) on the stored albedo, greenhouse quantity and the starlight
  at the orbit (a moon's is its planet's), floored at the cosmic microwave background; an airless body has G = 0.
- `ozone`: no oxygen, no ozone; with oxygen the column is at most the fit's median times the top of its spread (the abiotic branch's top
  below the fit), and at least the median times the bottom of the spread in a body that has not lost air.
- `aerosol-cloud`, `aerosol-haze`, `aerosol-dust`, `aerosol-airless`: the cloud depth and cover, the haze, the dust and the airless body's
  none follow `physics.aerosols`.
- `text-gases`: every gas the stored phrase names is in the composition table and the mixture, and the table's shares make one.
"""

import math

from planetgen.describe import build as describe_build
from planetgen.physics import aerosols, air_column, atmosphere, climate_draws, constants, gases, ozone, orbit_elements
from tests.audit.sweep import Rule, broken_copy
from tests.test_atmosphere_text_consistency import _text_names

REL = 1e-9
AU_KM = constants.AU_TO_KM


def _partials(body):
    return atmosphere.gas_partials_pa(body)


def _air(body):
    return sum(_partials(body).values()) > 0.0


def _gravity_m_s2(body):
    return body.gravity * constants.EARTH_GRAVITY


def _starlight_w_m2(entry):
    """The starlight at the body's orbit: a planet's at its own distance, a moon's at its planet's, on the mean light of the orbit."""
    host = entry.parent if entry.is_moon else entry.body
    eccentricity = entry.body.host_eccentricity if entry.is_moon else entry.body.eccentricity
    distance_au = orbit_elements.flux_distance_au(host.distance, eccentricity)
    return climate_draws.starlight_w_m2(entry.star.luminosity, distance_au * AU_KM)


def _close(value, expected, rel=REL):
    return value is not None and expected is not None and math.isclose(value, expected, rel_tol=rel, abs_tol=0.0)


def _mean_molar_mass(partials):
    total = sum(partials.values())
    return sum(pressure / total * gases.gas(formula).molar_mass for formula, pressure in partials.items())


# --- The rules ---------------------------------------------------------------------------------------------------------------


def _check_pressure_sum(entry):
    body, partials = entry.body, _partials(entry.body)
    if not partials:
        if body.atmospheric_pressure:
            return [f"no gases but a pressure of {body.atmospheric_pressure} Pa"]
        return []
    if not _close(sum(partials.values()), body.atmospheric_pressure):
        return [f"partials add up to {sum(partials.values())!r} Pa, stored pressure {body.atmospheric_pressure!r} Pa"]
    return []


def _break_pressure_sum(entry):
    if not _air(entry.body):
        return None
    return broken_copy(entry.body, atmospheric_pressure=entry.body.atmospheric_pressure * 1.01)


def _check_mean_molar_mass(entry):
    body = entry.body
    if not _air(body):
        return [] if body.atm_molar_density is None else [f"no air but a molar mass {body.atm_molar_density!r}"]
    expected = _mean_molar_mass(_partials(body))
    return [] if _close(body.atm_molar_density, expected) else [f"mean molar mass {body.atm_molar_density!r}, mixture gives {expected!r} kg/mol"]


def _break_mean_molar_mass(entry):
    return broken_copy(entry.body, atm_molar_density=entry.body.atm_molar_density * 1.01) if _air(entry.body) else None


def _check_air_density(entry):
    body = entry.body
    if not _air(body):
        return [] if body.atm_density is None else [f"no air but a density {body.atm_density!r}"]
    expected = body.atmospheric_pressure * _mean_molar_mass(_partials(body)) / (constants.R * body.surface_temperature)
    return [] if _close(body.atm_density, expected) else [f"density {body.atm_density!r}, P M / (R T) gives {expected!r} kg/m3"]


def _break_air_density(entry):
    return broken_copy(entry.body, atm_density=entry.body.atm_density * 1.01) if _air(entry.body) else None


def _check_scale_height(entry):
    body = entry.body
    if not _air(body):
        return [] if body.scale_height is None else [f"no air but a scale height {body.scale_height!r}"]
    expected = constants.R * body.surface_temperature / (_mean_molar_mass(_partials(body)) * _gravity_m_s2(body)) / 1000.0
    return [] if _close(body.scale_height, expected) else [f"scale height {body.scale_height!r}, R T / (M g) gives {expected!r} km"]


def _break_scale_height(entry):
    return broken_copy(entry.body, scale_height=entry.body.scale_height * 1.01) if _air(entry.body) else None


def _check_temperature(entry):
    body = entry.body
    if body.bond_albedo is None or body.greenhouse_factor is None:
        return ["no stored albedo or greenhouse quantity"]
    if not _air(body) and body.greenhouse_factor != 0.0:
        return [f"no air but a greenhouse quantity of {body.greenhouse_factor!r}"]
    sigma = constants.STEFAN_BOLTZMANN_CONSTANT
    flux = _starlight_w_m2(entry)
    expected = max(constants.COSMIC_BACKGROUND_TEMPERATURE_K,
                   ((1.0 - body.bond_albedo) * flux * (1.0 + body.greenhouse_factor) / (4.0 * sigma)) ** 0.25)
    return [] if _close(body.surface_temperature, expected) else [f"temperature {body.surface_temperature!r}, the albedo and greenhouse give {expected!r} K"]


def _break_temperature(entry):
    return broken_copy(entry.body, surface_temperature=entry.body.surface_temperature * 1.01)


def _o2_pa(body):
    return _partials(body).get("O2", 0.0)


def _ozone_bounds(entry):
    """`(low, high)` the stored column may take in DU for this body's oxygen and star, the low edge only for a body that kept its air."""
    body = entry.body
    pal = ozone.o2_level_pal(_o2_pa(body), body.gravity)
    teff = entry.star.temperature
    if pal <= 0.0:
        return 0.0, 0.0
    if pal < ozone.PAL_MIN:
        low, high = ozone.abiotic_column_range_du(teff)
        return low, high
    median = ozone.median_column_du(teff, pal)
    return median * ozone.STRUCTURE_RANGE[0], median * ozone.STRUCTURE_RANGE[1]


def _check_ozone(entry):
    column = getattr(entry.body, "ozone_column_du", None)
    if column is None:
        return []
    low, high = _ozone_bounds(entry)
    if column > high * (1.0 + REL) + 1e-12:
        return [f"ozone column {column!r} DU above the most the oxygen and star allow, {high!r}"]
    if column < low * (1.0 - REL) - 1e-12 and not _lost_air(entry):
        return [f"ozone column {column!r} DU below the least the oxygen and star allow, {low!r}"]
    return []


def _lost_air(entry):
    """Whether the body's ozone may have been cut by a loss event (a loss event that removes air); a sweep body has none."""
    return getattr(entry.body, "air_loss_events", 0) not in (0, None)


def _break_ozone(entry):
    column = getattr(entry.body, "ozone_column_du", None)
    if column is None:
        return None
    return broken_copy(entry.body, ozone_column_du=max(column * 1000.0, 1.0e4))


def _column(body):
    return air_column.column_of(body)


def _check_aerosol_airless(entry):
    body = entry.body
    if _air(body):
        return [f"{field} missing" for field in aerosols.BODY_FIELDS if getattr(body, field, None) is None]
    return [f"{field} = {getattr(body, field)!r} with no air" for field in aerosols.BODY_FIELDS if getattr(body, field, None) is not None]


def _break_aerosol_airless(entry):
    return broken_copy(entry.body, tau_haze=0.5) if not _air(entry.body) else broken_copy(entry.body, tau_haze=None)


def _check_aerosol_cloud(entry):
    body = entry.body
    column = _column(body)
    if column is None or body.tau_cloud is None:
        return []
    tau, cover = body.tau_cloud, body.cloud_cover
    most = aerosols.cloud_depth_limit(column, body.surface_temperature, aerosols.CLOUD_RADIUS_UM[0])
    if most <= 0.0:
        return [] if tau == 0.0 and cover == 0.0 else [f"no condensate but a cloud of depth {tau!r}, cover {cover!r}"]
    least = min(aerosols.CLOUD_TAU_FLOOR, aerosols.cloud_depth_limit(column, body.surface_temperature, aerosols.CLOUD_RADIUS_UM[1]))
    problems = []
    if not least * (1.0 - 1e-6) <= tau <= most * (1.0 + 1e-6):
        problems.append(f"cloud depth {tau!r} outside [{least!r}, {most!r}]")
    if not aerosols.CLOUD_COVER_RANGE[0] <= cover <= aerosols.CLOUD_COVER_RANGE[1]:
        problems.append(f"cloud cover {cover!r} outside {aerosols.CLOUD_COVER_RANGE}")
    return problems


def _break_aerosol_cloud(entry):
    column = _column(entry.body)
    if column is None or entry.body.tau_cloud is None:
        return None
    most = aerosols.cloud_depth_limit(column, entry.body.surface_temperature, aerosols.CLOUD_RADIUS_UM[0])
    return broken_copy(entry.body, tau_cloud=most * 10.0 + 1.0, cloud_cover=entry.body.cloud_cover or 0.5)


def _check_aerosol_haze(entry):
    body = entry.body
    column = _column(body)
    if column is None or body.tau_haze is None:
        return []
    if body.tau_haze < 0.0:
        return [f"negative haze depth {body.tau_haze!r}"]
    if aerosols.haze_gate(column) <= 0.0 and body.tau_haze != 0.0:
        return [f"haze depth {body.tau_haze!r} in air the methane to carbon dioxide ratio does not haze"]
    return []


def _break_aerosol_haze(entry):
    column = _column(entry.body)
    if column is None or entry.body.tau_haze is None:
        return None
    return broken_copy(entry.body, tau_haze=-1.0 if aerosols.haze_gate(column) > 0.0 else 0.5)


def _check_aerosol_dust(entry):
    body = entry.body
    if _column(body) is None or body.tau_dust is None:
        return []
    wet = body.body_type == "g" or bool(getattr(body, "ocean_fraction", None))
    if body.tau_dust == 0.0:
        return []
    low, high = aerosols.DUST_TAU_RANGE
    if wet:
        return [f"dust depth {body.tau_dust!r} over an ocean or on a giant"]
    return [] if low <= body.tau_dust <= high else [f"dust depth {body.tau_dust!r} outside {aerosols.DUST_TAU_RANGE}"]


def _break_aerosol_dust(entry):
    body = entry.body
    if _column(body) is None or body.tau_dust is None:
        return None
    return broken_copy(body, tau_dust=aerosols.DUST_TAU_RANGE[1] * 10.0)


def _check_text_gases(entry):
    body = entry.body
    if body.atmosphere in (None, "None") or not _air(body):
        return [] if body.atmosphere in (None, "None") and not _air(body) else [f"phrase {body.atmosphere!r} against gases {sorted(_partials(body))}"]
    if body.body_type == "g":
        return []
    block = describe_build._atmosphere_composition(body)
    table = {part["name"].lower(): part["share"] for part in block["parts"]}
    problems = []
    if not math.isclose(sum(table.values()), 1.0, abs_tol=1e-9):
        problems.append(f"table shares add up to {sum(table.values())!r}")
    mixture_names = {gases.gas(formula).name.lower() for formula in _partials(body)}
    if set(table) != mixture_names:
        problems.append(f"table gases {sorted(table)} against mixture {sorted(mixture_names)}")
    named = {name.lower() for name in _text_names(body.atmosphere)}
    if not named <= set(table):
        problems.append(f"phrase names {sorted(named - set(table))} that the table does not hold")
    return problems


def _break_text_gases(entry):
    body = entry.body
    if body.body_type == "g" or body.atmosphere in (None, "None") or not _air(body):
        return None
    return broken_copy(body, atmosphere="a thin mix of unobtainium")


RULES = (
    Rule("pressure-sum", _check_pressure_sum, _break_pressure_sum),
    Rule("mean-molar-mass", _check_mean_molar_mass, _break_mean_molar_mass),
    Rule("air-density", _check_air_density, _break_air_density),
    Rule("scale-height", _check_scale_height, _break_scale_height),
    Rule("temperature", _check_temperature, _break_temperature),
    Rule("ozone", _check_ozone, _break_ozone),
    Rule("aerosol-airless", _check_aerosol_airless, _break_aerosol_airless),
    Rule("aerosol-cloud", _check_aerosol_cloud, _break_aerosol_cloud),
    Rule("aerosol-haze", _check_aerosol_haze, _break_aerosol_haze),
    Rule("aerosol-dust", _check_aerosol_dust, _break_aerosol_dust),
    Rule("text-gases", _check_text_gases, _break_text_gases),
)
