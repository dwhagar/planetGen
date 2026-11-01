# planetgen/physics/giant_radius.py

"""
The radius of a giant planet from its mass, heavy elements, irradiation and age (GEN.460), and the heavy-element mass of
a giant [R0710].

**Radius.** The tables of [R1429] (`giant_radius_tables`, Tables 2 to 4 of ApJ 659:1661) give a
cool-start radius for masses from 17 to 3,594 Earth masses, five distances from a one-solar-luminosity star, five core
masses and three ages. `radius_rj` interpolates them: linearly in the core mass, and in the logarithms of mass, distance
(the flux, `d = 1/sqrt(S)` for a flux `S` in Earth units) and age. Outside the grid the nearest edge is held, in flux
and age, where a giant's radius changes slowly; in mass the radius is held below 17 Earth masses (where `in_range` is
false and the caller uses another model) and continued with the last two columns up to `MASS_REACH` times the heaviest.
The table's core is every heavy element the planet has outside its envelope: `table_core_mass_earth` takes the envelope's
own heavy elements (its metallicity) off the planet's `M_Z`. The tables put the whole core in the middle, so a planet
whose heavy elements are spread through its interior (a diluted core) is larger than the table says: Jupiter, with 42
Earth masses of heavy elements, comes out at 0.95 Jupiter radii. A planet that is mostly heavy elements (an ice giant)
is outside the tables, and `covers` says so.

**Inflation.** The tables do not inflate. Above a threshold flux a giant of 0.37 Jupiter masses or more is larger than
a cool one, and `inflation_rj` adds the population relation of [R4029] (A&A 612:A74, Table 1),
`Delta R = A (log10 F - log10 F_s)` Jupiter radii for a flux `F` in W/m^2 above `F_s`, in three mass bins. It is an
empirical fit to 286 transiting giants, not an interior model: the physical efficiency of [R4351]
(ApJ 155:214, 2.4 percent of the flux at 1,500 K, falling to 0.2 percent at 2,500 K) needs a thermal evolution code the
project does not have, so it stays a stand-in (marked in the gap plan). The relation is the median; the fraction of
hot giants that are not inflated is not drawn.

**Heavy elements.** `heavy_mass_earth` draws `M_Z = M_core + f_Z (M - M_core)` ([R0710],
equation 6; `giant-heavy-element-law.md`) with the paper's multiplicative scatter, and `M_Z = M` below `M_core`, kept
below `MAX_HEAVY_FRACTION` of the mass so every giant keeps an envelope.
"""

import bisect
import math

from planetgen.physics import constants, giant_radius_tables as tables, unit_registry as units
from planetgen.util import draw

_kg_to_mearth = units.converter("kg", "Mearth")
_kg_to_mjup = units.converter("kg", "Mjup")
_mearth_to_kg = units.converter("Mearth", "kg")
_mjup_to_kg = units.converter("Mjup", "kg")

MASS_REACH = 1.35
"""float: How far past the heaviest tabulated mass (3,594 Earth masses) the last two columns are continued, as a
factor; a brown dwarf at the deuterium limit is 1.2 times it."""

MASSES_EARTH = tuple(_kg_to_mearth(_mjup_to_kg(m)) for m in tables.MASSES_MJ)
"""tuple: The table's masses in Earth masses."""

INFLATION_BINS = (
    # (lowest mass in Jupiter masses, A in Jupiter radii per dex, log10 of the threshold flux in W/m^2)
    (0.37, 0.70, 5.52),
    (0.98, 0.52, 5.82),
    (2.50, 0.22, 5.20),
)
"""tuple: [R4029: Table 1]: the mass bins from 0.37 Jupiter masses, with `A` and `log10 F_s`. Below 0.37
Jupiter masses the paper's `A` is -0.33 (+0.51, -0.91), no inflation within its error, so none is applied."""

CHACHAN_CORE_MASS_EARTH = 14.7
"""float: `M_core` of [R0710] (+1.8, -1.6 Earth masses)."""

CHACHAN_ENVELOPE_HEAVY_SHARE = 0.09
"""float: `f_Z` of [R0710] (+-0.01): the heavy-element share of the mass accreted beyond `M_core`."""

CHACHAN_SCATTER = 0.66
"""float: `sigma_mult` of [R0710] (+-0.08): the spread of `M_Z` as a fraction of its predicted value."""

MAX_HEAVY_FRACTION = 0.9
"""float: The most of a giant's mass that is heavy elements; the rest is hydrogen and helium."""

MIN_HEAVY_FRACTION = 0.01
"""float: The least heavy-element share a draw may have, about the Sun's (0.014), since an envelope at least as rich as
its star's is the lower bound of the scatter."""

DRAW_ATTEMPTS = 64
"""int: How many scatter draws are tried before the nearest bound is taken."""


class GiantRadiusError(ValueError):
    """A mass, flux or age the tables cannot be read at."""


def table_core_mass_earth(mass_earth, heavy_mass_earth, envelope_z=0.0):
    """The core mass the tables are read at: the heavy elements less the share of them in the envelope, which has the
    heavy-element mass fraction `envelope_z` of the mass outside the core."""
    return max(0.0, (heavy_mass_earth - envelope_z * mass_earth) / (1.0 - envelope_z))


def covers(mass_earth, core_mass_earth):
    """Whether the tables hold a core as heavy as `core_mass_earth` at `mass_earth`: the heaviest core they have a value
    for at the tabulated mass just below it."""
    below = max(i for i, m in enumerate(MASSES_EARTH) if m <= max(mass_earth, MASSES_EARTH[0]) * (1.0 + 1e-9))
    column = tables.RADII_RJ[-1][-1]  # the 4.5 Gy, 9.5 AU table has a value wherever the paper has one
    heaviest = max(c for c, radii in zip(tables.CORE_MASSES_EARTH, column) if radii[below] is not None)
    return core_mass_earth <= heaviest


def in_range(mass_earth):
    """Whether the tables (with `MASS_REACH`) cover `mass_earth`, in Earth masses."""
    return MASSES_EARTH[0] * (1.0 - 1e-9) <= mass_earth <= MASSES_EARTH[-1] * MASS_REACH


def _bracket(grid, value):
    """`(index, share)` of `value` between `grid[index]` and `grid[index + 1]`, held at the ends."""
    if value <= grid[0]:
        return 0, 0.0
    if value >= grid[-1]:
        return len(grid) - 2, 1.0
    index = bisect.bisect_right(grid, value) - 1
    return index, (value - grid[index]) / (grid[index + 1] - grid[index])


def _by_core(radii, core_earth):
    """The radius of one `[core]` column at `core_earth`, interpolated over the cores the table has a value for and held
    at the lightest and heaviest of those."""
    known = [(c, r) for c, r in zip(tables.CORE_MASSES_EARTH, radii) if r is not None]
    if len(known) == 1:
        return known[0][1]
    index, share = _bracket([c for c, _ in known], core_earth)
    return known[index][1] + share * (known[index + 1][1] - known[index][1])


def _cell(age_index, distance_index, mass_index, core_earth):
    column = [tables.RADII_RJ[age_index][distance_index][c][mass_index] for c in range(len(tables.CORE_MASSES_EARTH))]
    return _by_core(column, core_earth)


def _at_mass(age_index, distance_index, mass_earth, core_earth):
    """The radius at one age and distance node: log-linear in mass, continued past the last column."""
    logs = [math.log(m) for m in MASSES_EARTH]
    x = math.log(max(mass_earth, MASSES_EARTH[0]))
    if x > logs[-1]:
        low = _cell(age_index, distance_index, len(logs) - 2, core_earth)
        high = _cell(age_index, distance_index, len(logs) - 1, core_earth)
        return high + (x - logs[-1]) / (logs[-1] - logs[-2]) * (high - low)
    index, share = _bracket(logs, x)
    low = _cell(age_index, distance_index, index, core_earth)
    high = _cell(age_index, distance_index, index + 1, core_earth)
    return low + share * (high - low)


def radius_rj(mass_earth, core_mass_earth, flux_earth, age_gy):
    """
    The radius in Jupiter radii (equatorial, 1 bar) of a giant of `mass_earth` Earth masses with a core of
    `core_mass_earth` Earth masses (`table_core_mass_earth`), irradiated at `flux_earth` times the Earth's flux, at
    `age_gy` Gy. Without inflation. A core heavier than the table has one for this mass is held at the heaviest it has
    (`covers` is false then, and the radius is too large for a planet that rich in heavy elements).

    Raises:
        GiantRadiusError: The mass is outside `in_range`, or the flux or age is not positive.
    """
    if not in_range(mass_earth):
        raise GiantRadiusError(f"a {mass_earth} Earth mass giant is outside the tables ({MASSES_EARTH[0]:.0f} to "
                               f"{MASSES_EARTH[-1] * MASS_REACH:.0f})")
    if flux_earth <= 0.0 or age_gy <= 0.0:
        raise GiantRadiusError("the flux and the age must be positive")
    log_distance = -0.5 * math.log(flux_earth)
    distances = [math.log(d) for d in tables.DISTANCES_AU]
    ages = [math.log(a) for a in tables.AGES_GY]
    di, dshare = _bracket(distances, log_distance)
    ai, ashare = _bracket(ages, math.log(age_gy))
    total = 0.0
    for age_node, aw in ((ai, 1.0 - ashare), (ai + 1, ashare)):
        for distance_node, dw in ((di, 1.0 - dshare), (di + 1, dshare)):
            if aw and dw:
                total += aw * dw * _at_mass(age_node, distance_node, mass_earth, core_mass_earth)
    return total


def inflation_rj(mass_earth, flux_w_m2):
    """The radius a hot giant gains, in Jupiter radii, from the population relation of [R4029] (see the
    module text); 0 below 0.37 Jupiter masses and below the mass bin's threshold flux."""
    mass_mj = _kg_to_mjup(_mearth_to_kg(mass_earth))
    bin_ = None
    for row in INFLATION_BINS:
        if mass_mj >= row[0]:
            bin_ = row
    if bin_ is None or flux_w_m2 <= 0.0:
        return 0.0
    _, amplitude, log_threshold = bin_
    return max(0.0, amplitude * (math.log10(flux_w_m2) - log_threshold))


def flux_w_m2(flux_earth):
    """A flux in Earth units as watts per square metre (the Earth's is the Sun's luminosity over 4 pi AU^2)."""
    return flux_earth * constants.SOLAR_LUMINOSITY / (4.0 * math.pi * constants.AU_M ** 2)


def radius_km(mass_earth, core_mass_earth, flux_earth, age_gy):
    """`radius_rj` plus `inflation_rj`, in km."""
    radius = radius_rj(mass_earth, core_mass_earth, flux_earth, age_gy) + inflation_rj(mass_earth, flux_w_m2(flux_earth))
    return radius * constants.JUPITER_RADIUS_KM


def predicted_heavy_mass_earth(mass_earth):
    """[R0710]'s heavy-element mass of a giant of `mass_earth`, before scatter: `M_Z = M_core + f_Z (M -
    M_core)`, and all of the mass below `M_core`."""
    if mass_earth <= CHACHAN_CORE_MASS_EARTH:
        return mass_earth
    return CHACHAN_CORE_MASS_EARTH + CHACHAN_ENVELOPE_HEAVY_SHARE * (mass_earth - CHACHAN_CORE_MASS_EARTH)


def heavy_mass_earth(mass_earth, rng=draw):
    """
    A giant's heavy-element mass in Earth masses: the predicted mass (`predicted_heavy_mass_earth`) times `1 + g`, with
    `g` normal of standard deviation `CHACHAN_SCATTER`, redrawn until the result lies between `MIN_HEAVY_FRACTION` and
    `MAX_HEAVY_FRACTION` of the mass (the nearest bound after `DRAW_ATTEMPTS`). The paper gives the spread as a fraction
    of the prediction, not its shape; a normal one is the stand-in.
    """
    predicted = predicted_heavy_mass_earth(mass_earth)
    low, high = MIN_HEAVY_FRACTION * mass_earth, MAX_HEAVY_FRACTION * mass_earth
    value = predicted
    for _ in range(DRAW_ATTEMPTS):
        value = predicted * (1.0 + CHACHAN_SCATTER * rng.gauss(0.0, 1.0))
        if low <= value <= high:
            return value
        draw.redrawn("giant.heavy_mass")
    return min(max(value, low), high)
