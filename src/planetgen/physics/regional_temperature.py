# planetgen/physics/regional_temperature.py

"""
Regional temperatures and the day and night swing (GEN.656; docs/design/habitability/survival-equipment-fields.md,
section 5.5, Boss's decision C of 2026-10-23).

The PHI-6 score and the required survival gear read one average surface temperature. This module spreads that mean
over the world so the detail text can say where it is colder or hotter than the mean, and by how much the day and the
night differ. Everything is a pure function of stored columns (no draws, no epoch):

* **Latitude bands** (a world that turns on its own): Polar (60 degrees and over, 13 percent of the surface),
  Northern (20 to 60 N, 26 percent), Equatorial (20 S to 20 N, 34 percent), Southern (20 to 60 S, 26 percent). Each
  takes the annual-mean insolation of its band on a circular orbit for the stored axial tilt (`band_insolation`),
  turned into a temperature by the fourth-root law and scaled so the area-weighted mean of the bands is the world's
  mean. The air then carries heat poleward (`transport_share`): the share `x / (x + TRANSPORT_X0)` of the contrast
  is removed, `x` being the air column (pressure over gravity) in Earth's.
* **Tidally locked worlds** (the rotation period is the orbital period): the day side, the terminator ring and the
  night side take the place of the bands, and the same transport share pulls them toward the mean.
* **The swing**: a damped-oscillator estimate. Without damping the equatorial noon is `sqrt(2)` times the mean and
  the night is cold; the heat capacity of the air column (`c_p P / g`) plus the skin of the ground
  (`I sqrt(period / pi)`) against the radiative loss `4 sigma T^3` gives a time constant `tau`, and the swing is cut
  by `sqrt(1 + (omega tau)^2)` and by the share of the contrast the air carries away.

The orbits here are circular, so Northern and Southern are equal; the columns stay separate for the eccentric orbit
the model will gain.
"""

import math

from planetgen.physics import constants

FIELDS = (
    "temp_polar_k", "temp_north_k", "temp_equator_k", "temp_south_k", "temp_swing_k",
    "temp_day_side_k", "temp_terminator_k", "temp_night_side_k",
)
"""tuple: The stored columns on `planets` and `moons`: the four latitude bands and the swing for a world that turns on
its own, the three zones for a tidally locked one (the other group is `None`), all in kelvin."""

BANDS = (("polar", 60.0, 90.0), ("north", 20.0, 60.0), ("equator", -20.0, 20.0), ("south", -60.0, -20.0))
"""tuple: `(name, from latitude, to latitude)` of the bands, degrees. Polar is the cap above 60 degrees at either pole
(the same on a circular orbit)."""

TRANSPORT_X0 = 20.0
"""float: The air column, in Earth's, at which the air removes half of the pole to equator contrast (DESIGN). Fitted
so that Earth's polar band sits about 35 K under its mean (the observed annual polar mean is 35 to 40 K under), which
makes the share 0.05 for Earth, 0.3 for Titan and 0.8 for Venus."""

LOCKED_TRANSPORT_X0 = 0.1
"""float: The same for the day to night contrast of a tidally locked world (DESIGN): a day side with no turning to
spread the heat needs far less air to share it, so a bar of air removes 90 percent of the contrast (the order the
general circulation models of locked Earth-like worlds give) and a Mars-thin one 14 percent."""

EARTH_COLUMN_PA_PER_MS2 = 101_325.0 / constants.EARTH_GRAVITY
"""float: Earth's air column, pressure over gravity, Pa s2/m: the unit of `x`."""

AIR_CP_J_KG_K = 1000.0
"""float: The heat capacity of air at constant pressure, J/(kg K) (about 1004 for Earth's, 850 for carbon dioxide,
14 000 for hydrogen; one value for all, a DESIGN stand-in until the mixture's own is wired)."""

GROUND_INERTIA = 300.0
"""float: The thermal inertia of the ground, J/(m2 K s^0.5) (DESIGN: lunar regolith 55, Mars about 300, bare rock
2000 to 2500; a stand-in until the surface composition sets it)."""

NIGHT_FLOOR_K = 20.0
"""float: The coldest a night side drops, kelvin (the heat from the ground; DESIGN stand-in)."""

LOCK_TOLERANCE = 1e-3
"""float: How close the rotation period is to the orbital period to count as locked (the generator sets them equal)."""

TERMINATOR_RADIATIVE = 0.5
"""float: The terminator ring, as a fraction of the mean, for a locked world without air (DESIGN)."""

DAY_SIDE_RADIATIVE = 1.131
"""float: The day side, as a multiple of the mean, for a locked world without air: the hemisphere average of the
fourth root of the cosine of the angle from the sub-stellar point is 4/5, times the sub-stellar `4^(1/4)`."""

SUBSTELLAR_RATIO = math.sqrt(2.0)
"""float: The noon equator over the mean for a world that cools completely at night (`4^(1/4)`)."""

_LATITUDE_STEPS = 24
_SEASON_STEPS = 72


def _daily_insolation(latitude, declination):
    """The day's mean insolation at `latitude` for the Sun at `declination`, in units of the solar constant (radians)."""
    x = -math.tan(latitude) * math.tan(declination)
    if x >= 1.0:
        return 0.0
    hour = math.pi if x <= -1.0 else math.acos(x)
    return (hour * math.sin(latitude) * math.sin(declination)
            + math.cos(latitude) * math.cos(declination) * math.sin(hour)) / math.pi


def _annual_insolation(latitude, tilt):
    """The year's mean at `latitude` on a circular orbit with axial tilt `tilt` (radians), units of the solar constant."""
    total = 0.0
    for step in range(_SEASON_STEPS):
        longitude = (step + 0.5) / _SEASON_STEPS * 2.0 * math.pi
        total += _daily_insolation(latitude, math.asin(math.sin(tilt) * math.sin(longitude)))
    return total / _SEASON_STEPS


def band_insolation(tilt_deg):
    """
    `{band: annual-mean insolation over the planet's mean}` for the four bands of `BANDS` (area-weighted over the
    band's latitudes, so a sphere averages 1). The south pole's cap and the north pole's are the same by symmetry of
    the circular orbit, so the Polar band is the north cap.
    """
    tilt = math.radians(min(max(tilt_deg or 0.0, 0.0), 180.0))
    out = {}
    for name, low, high in BANDS:
        weight = total = 0.0
        for step in range(_LATITUDE_STEPS):
            latitude = math.radians(low + (step + 0.5) / _LATITUDE_STEPS * (high - low))
            weight += math.cos(latitude)
            total += math.cos(latitude) * _annual_insolation(abs(latitude), tilt)
        out[name] = 4.0 * total / weight
    return out


def band_areas():
    """`{band: share of the surface}`: the Polar cap counts both poles (13 percent)."""
    shares = {}
    for name, low, high in BANDS:
        share = (math.sin(math.radians(high)) - math.sin(math.radians(low))) / 2.0
        shares[name] = share * (2.0 if name == "polar" else 1.0)
    return shares


def column_in_earths(pressure_pa, gravity_ms2):
    """The air column (pressure over gravity) in Earth's; 0 for no air."""
    if not pressure_pa or not gravity_ms2 or pressure_pa <= 0.0 or gravity_ms2 <= 0.0:
        return 0.0
    return pressure_pa / gravity_ms2 / EARTH_COLUMN_PA_PER_MS2


def transport_share(column, locked=False):
    """The share of the contrast the air removes, 0 (no air) toward 1 (a deep atmosphere), from the column in Earth's."""
    x0 = LOCKED_TRANSPORT_X0 if locked else TRANSPORT_X0
    return column / (column + x0) if column > 0.0 else 0.0


def is_locked(rotation_period_hours, period_years):
    """Whether the body turns once an orbit (the generator's lock): its rotation period is its orbital period."""
    if not rotation_period_hours or not period_years:
        return False
    orbital_hours = period_years * constants.SECONDS_PER_YEAR / 3600.0
    return abs(rotation_period_hours - orbital_hours) <= LOCK_TOLERANCE * orbital_hours


def swing_k(mean_k, pressure_pa, gravity_ms2, rotation_period_hours):
    """
    The day and night swing, peak to peak at the equator, kelvin: `sqrt(2) T` cut by the damping of the air column and
    the ground skin (see the module text). A fast turn with a thick air column gives a few kelvin, a slow turn on bare
    ground the full `sqrt(2) T`.
    """
    omega = 2.0 * math.pi / (rotation_period_hours * 3600.0)
    column = (pressure_pa / gravity_ms2) if pressure_pa and gravity_ms2 and pressure_pa > 0.0 else 0.0
    capacity = AIR_CP_J_KG_K * column + GROUND_INERTIA * math.sqrt(1.0 / omega * 2.0)
    time_constant = capacity / (4.0 * constants.STEFAN_BOLTZMANN_CONSTANT * mean_k ** 3)
    return SUBSTELLAR_RATIO * mean_k / math.sqrt(1.0 + (omega * time_constant) ** 2)


def compute(values):
    """
    The `FIELDS` of a rocky body from its column-named `values` (`habitability_world.values_from_body`, or a row):
    `surface_temperature_k`, `atmospheric_pressure_pa`, `gravity_g`, `rotation_period_hours`, `axial_tilt_deg` and
    `period_years`. Everything is `None` for a body without a temperature or a rotation period.
    """
    out = dict.fromkeys(FIELDS)
    mean = values.get("surface_temperature_k")
    rotation = values.get("rotation_period_hours")
    if not mean or mean <= 0.0 or not rotation or rotation <= 0.0:
        return out
    pressure = values.get("atmospheric_pressure_pa") or 0.0
    gravity = (values.get("gravity_g") or 0.0) * constants.EARTH_GRAVITY
    column = column_in_earths(pressure, gravity)
    locked = is_locked(rotation, values.get("period_years"))
    share = transport_share(column, locked)
    if locked:
        zones = {"day_side": DAY_SIDE_RADIATIVE, "terminator": TERMINATOR_RADIATIVE, "night_side": 0.0}
        for name, radiative in zones.items():
            temperature = mean + (radiative * mean - mean) * (1.0 - share)
            out[f"temp_{name}_k"] = max(temperature, NIGHT_FLOOR_K)
        return out
    insolation = band_insolation(values.get("axial_tilt_deg"))
    areas = band_areas()
    ratios = {name: insolation[name] ** 0.25 for name in insolation}
    scale = sum(areas[name] * ratios[name] for name in ratios) / sum(areas.values())
    for name in ratios:
        radiative = ratios[name] / scale
        out[f"temp_{name}_k"] = mean + (radiative * mean - mean) * (1.0 - share)
    out["temp_swing_k"] = swing_k(mean, pressure, gravity, rotation) * (1.0 - share)
    return out


def summary(values):
    """The regions object the API serves, kelvin by short name (`polar`, `north`, `equator`, `south`, `swing`, `day_side`,
    `terminator`, `night_side`; only the ones the body has), or `None` for a body that has none."""
    out = {field[len("temp_"):-len("_k")]: values[field] for field in FIELDS if values.get(field) is not None}
    return out or None
