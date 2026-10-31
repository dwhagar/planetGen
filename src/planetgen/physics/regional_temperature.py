# planetgen/physics/regional_temperature.py

"""
Regional temperatures and the day and night swing (GEN.656, GEN.1108; docs/design/habitability/survival-equipment-fields.md,
section 5.5, Boss's decision C of 2026-10-23).

The PHI-6 score and the required survival gear read one average surface temperature. This module spreads that mean
over the world so the detail text can say where it is colder or hotter than the mean, and by how much the day and the
night differ. Everything is a pure function of stored columns (no draws, no epoch). Since GEN.1108 the places are the seven
bands of `physics/latitude_bands.py` (the equator band, three on each side of it, the poles inside the outermost): for each
`band1` (north pole) to `band7` (south pole) the stored columns are `bandN_mean_k`, `bandN_swing_k`, `bandN_day_hours`,
`bandN_night_hours`, `bandN_longest_day_days`, `bandN_longest_night_days` and `bandN_state` (`column`), and the world's mean
swing is `temp_mean_swing_k`. A tidally locked world has no latitude day: its seven bands are the rings from the point under
its star to the point away from it (`band1` the day centre, `band7` the night centre), with no hours.

* **Temperatures**: each band takes the annual-mean insolation of its band on a circular orbit for the stored axial tilt
  (`band_insolation`), turned into a temperature by the fourth-root law and scaled so the area-weighted mean of the bands is the
  world's mean. The air then carries heat poleward (`transport_share`): the share `x / (x + TRANSPORT_X0)` of the contrast
  is removed, `x` being the air column (pressure over gravity) in Earth's.
* **The swing**: a damped-oscillator estimate. Without damping the equatorial noon is `sqrt(2)` times the mean and
  the night is cold; the heat capacity of the air column (`c_p P / g`) plus the skin of the ground
  (`I sqrt(period / pi)`) against the radiative loss `4 sigma T^3` gives a time constant `tau`, and the swing is cut
  by `sqrt(1 + (omega tau)^2)` and by the share of the contrast the air carries away; a band's is the equator's times the
  cosine of its latitude.

The orbits here are circular, so a north band and its south twin are equal; the columns stay separate for the eccentric orbit
the model will gain.
"""

import math

import numpy as np

from planetgen.physics import constants
from planetgen.physics import unit_registry as units

_yr_to_s = units.converter("yr", "s")
_s_to_h = units.converter("s", "h")

BAND_COUNT = 7
"""int: How many bands a body's climate is stored in (`latitude_bands`)."""

BAND_FIELDS = ("mean_k", "swing_k", "day_hours", "night_hours", "longest_day_days", "longest_night_days", "state")
"""tuple: What is stored for each band: the mean temperature (K), the day and night swing, peak to peak (K), the hours of day and of night in a solar day, the longest
polar day and night of the year (days) at the poleward edge, and the state word."""


def column(number, field):
    """The stored column of band `number` (1, the north pole, to 7, the south pole; for a locked world 1 the day centre to 7 the night centre) and `field` of `BAND_FIELDS`."""
    return f"band{number}_{field}"


BAND_COLUMNS = tuple(column(number, field) for number in range(1, BAND_COUNT + 1) for field in BAND_FIELDS)
"""tuple: Every band column, band by band."""

FIELDS = (*BAND_COLUMNS, "temp_mean_swing_k")
"""tuple: The stored columns on `planets` and `moons` (kelvin, hours, days and a state word): the seven bands, and the mean swing of the world (the area-weighted mean of
the bands'). All `None` for a body with no temperature or rotation period."""

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

_SEASON_STEPS = 72


def _daily_insolation(latitude, declination):
    """The day's mean insolation at `latitude` (array, radians) for the Sun at `declination` (array), in units of the
    solar constant: 0 in the polar night, the whole day's arc in the polar day."""
    x = -np.tan(latitude) * np.tan(declination)
    hour = np.arccos(np.clip(x, -1.0, 1.0))                  # pi for x <= -1 (the sun never sets), 0 for x >= 1
    return (hour * np.sin(latitude) * np.sin(declination)
            + np.cos(latitude) * np.cos(declination) * np.sin(hour)) / math.pi


def _annual_insolation(latitude, tilt):
    """The year's mean at `latitude` (array, radians) on a circular orbit with axial tilt `tilt` (radians), units of the
    solar constant: the mean of the day's over `_SEASON_STEPS` points of the year."""
    longitude = (np.arange(_SEASON_STEPS) + 0.5) / _SEASON_STEPS * 2.0 * math.pi
    declination = np.arcsin(math.sin(tilt) * np.sin(longitude))
    return _daily_insolation(np.asarray(latitude)[:, None], declination[None, :]).mean(axis=1)


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
    orbital_hours = _s_to_h(_yr_to_s(period_years))
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


def daily_range(value_k, swing_k):
    """
    `(low, high)` kelvin: the coldest and hottest hour of a place whose stored temperature is `value_k` and whose day and
    night swing, peak to peak, is `swing_k`, the value exactly in the middle (GEN.795, GEN.796). The swing is a
    difference, so half of it is added to and taken from the temperature; it is never converted as a temperature.
    `None` when either is missing or the swing is 0.
    """
    if value_k is None or not swing_k:
        return None
    half = swing_k / 2.0
    return value_k - half, value_k + half


def compute(values):
    """
    The `FIELDS` of a rocky body from its column-named `values` (`habitability_world.values_from_body`, or a row):
    `surface_temperature_k`, `atmospheric_pressure_pa`, `gravity_g`, `rotation_period_hours`, `axial_tilt_deg`, `period_years` and the
    mean share of the time the star is up, `day_share_mean`, from which the angular radius of the star is worked out. Everything is
    `None` for a body without a temperature or a rotation period.
    """
    from planetgen.physics import latitude_bands       # here: it reads this module's laws
    out = dict.fromkeys(FIELDS)
    bands = latitude_bands.compute(values, latitude_bands.alpha_deg_of(values))
    if not bands:
        return out
    for number, band in enumerate(bands, 1):
        for field in BAND_FIELDS:
            out[column(number, field)] = getattr(band, field)
    out["temp_mean_swing_k"] = sum(band.area_share * (band.swing_k or 0.0) for band in bands) / sum(band.area_share for band in bands)
    return out


def summary(values):
    """
    The regions object the API serves (`None` for a body that has none): `bands`, the seven bands from the north pole to the south pole (the day centre to the night centre
    for a world that turns once an orbit), each `{"name", "mean_k", "swing_k", "day_hours", "night_hours", "longest_day_days", "longest_night_days", "state"}` with the
    fields the body has, and `mean_swing_k`. Names are `latitude_bands.NAMES`, or `RING_NAMES` for a locked world.
    """
    from planetgen.physics import latitude_bands
    bands = []
    for number in range(1, BAND_COUNT + 1):
        entry = {field: values.get(column(number, field)) for field in BAND_FIELDS}
        if entry["mean_k"] is None:
            return None
        state = entry["state"]
        names = latitude_bands.RING_NAMES if state in latitude_bands.RING_STATES else latitude_bands.NAMES
        bands.append({"name": names[number - 1], **{key: value for key, value in entry.items() if value is not None}})
    return {"bands": bands, "mean_swing_k": values.get("temp_mean_swing_k")}
