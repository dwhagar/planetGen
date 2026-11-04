# planetgen/galaxy/disk_history.py

"""
The Thin Disk's Star Formation History (GEN.469)
================================================

How the thin disk's stars are spread over their ages, worked out from the disk's own gas instead of typed in
(`docs/design/galaxy/disk-star-formation-history.md`). The disk's three populations (`young`, `intermediate`, `old`) share its
stars in proportion to how much star formation fell in their age windows; a flat history gave each window's length over
the disk's age, and this module gives what the gas gives.

Each ring of the disk (a thousandth of its stars or so) follows gas continuity:

    d(gas)/dt = inflow(t) - (1 - R) * rate(gas),    rate = A * gas^n

1. **The depletion** is the Kennicutt-Schmidt law ([R2321], `tuning.STAR_FORMATION_LAW_COEFFICIENT` and
   `STAR_FORMATION_LAW_EXPONENT`): the star formation per area goes as the gas per area to the 1.4, steepened where the
   gas thins by the factor `(gas / (gas + saturation)) ^ (m - n)` of GEN.450 (`tuning.STAR_FORMATION_SATURATION_GAS_MSUN_PC2`,
   `STAR_FORMATION_LOW_DENSITY_EXPONENT`; GEN.999). Without it a thin outer ring spends its gas as fast as a dense one and a
   Milky Way mass disc ended with a third of the gas it has.
2. **The return** `R` is the share of a stellar generation's mass that goes back into the gas (`return_fraction`),
   worked out from the generator's own mass function with the remnant masses of [R2798].
3. **The inflow** is an exponential infall whose time scale grows outward (inside-out formation): `tau(R) = 1.033 R - 1.267` Gy of the
   Milky Way's chemical evolution model ([R5095], section 3), stopped at 2 Gy inside and 8 Gy outside, the shortest and longest of
   the later model of [R5094] (`infall_time_gy`, `tuning.DISK_INFALL_*`). The ring at `x` scale lengths has the time scale of the
   Milky Way's ring at `x` of its own 2.6 kpc: no source gives how the infall scales with a galaxy's size (GEN.725 row GI1). The disc's mass enters through
   its circular speed `V` (the cube root of its mass): the inflow is also multiplied by `exp(0.4 t (1 - V / 220))` per Gy ([R5096] equation 21), so a lighter
   disc is fed more and more slowly and forms its stars late, a heavier one the other way (GEN.1246).
   Only its shape in time matters: each ring's amplitude is solved so that the ring ends with the stars the galaxy's density
   gives it. (Until GEN.1230 the inflow followed the stars the galaxy's halo builds, which for the lighter halos of GEN.1154 kept
   the disk forming 2.7 times the stars per year the main sequence gives.)
4. **The stars.** The disk's stellar mass (its share of the galaxy's, which the mass chain drew from the seed) spread as
   `exp(-R / L)` sets each ring's stellar surface density, which fixes how fast its gas is used (`rate / gas` is
   `A * gas^(n-1)`): a dense ring follows the inflow closely, a thin one lags it.

What is left out, as in any one-zone model: gas moving between rings, outflows, and the stars a halo gains by mergers
(they are in the relation's stellar mass, not in the thin disk). The ring model is run for the disk-wide time shape
only; the radial pattern of recent star formation stays `galaxy.star_formation`'s (GEN.134), because the rings alone
make the outer disk younger, against the observed profile.
"""

import dataclasses
import functools
import math

import numpy as np

from planetgen import tuning
from planetgen.galaxy import mass_chain
from planetgen.galaxy.age_law import AgeLaw
from planetgen.physics import imf, mass_boundaries, remnant_mass

_BISECTION_STEPS = 40
_DOUBLINGS = 80


@dataclasses.dataclass(frozen=True)
class DiskHistory:
    """
    The rings' histories.

    Attributes:
        radius_scale_lengths (numpy.ndarray): Each ring's middle, in disk scale lengths.
        age_gy (numpy.ndarray): The time before now of each step, decreasing from the disk's age to 0.
        rate (numpy.ndarray): Star formation per area, rings by steps, solar masses per Gy per pc2 (births).
        gas_now (numpy.ndarray): Each ring's gas per area at the end, solar masses per pc2.
        stars_now (numpy.ndarray): Each ring's stellar surface density it was solved to reach, solar masses per pc2.
    """
    radius_scale_lengths: np.ndarray
    age_gy: np.ndarray
    rate: np.ndarray
    gas_now: np.ndarray
    stars_now: np.ndarray

    def share_between(self, young_edge_gy, old_edge_gy):
        """
        The share of the disk's stars (mass formed, rings weighted by their area) born between `young_edge_gy` and
        `old_edge_gy` before now.
        """
        span = float(self.age_gy[0])
        edges = np.clip(np.array([young_edge_gy, old_edge_gy], dtype=float), 0.0, span)
        total_rate = self._total_rate()
        step = span / (len(self.age_gy) - 1)
        cumulative = np.concatenate(([0.0], np.cumsum(0.5 * (total_rate[1:] + total_rate[:-1]) * step)))
        since_start = span - edges                               # time since the disk began, at each edge
        at = np.interp(since_start, span - self.age_gy, cumulative)
        return float((at[0] - at[1]) / cumulative[-1])

    def age_law(self, young_edge_gy, old_edge_gy, points=None):
        """
        The spread of ages inside the window from `young_edge_gy` to `old_edge_gy` before now (GEN.546): the disk's star
        formation (rings weighted by their area) at `points` ages spread evenly across it, as an `AgeLaw`. The rates are
        rounded to six figures so that two histories that agree to that give one law (the sampling tables are cached on
        it).
        """
        points = tuning.AGE_LAW_POINTS if points is None else points
        ages = np.linspace(young_edge_gy, old_edge_gy, points)
        # `age_gy` falls from the disk's age to 0, so the rate is read on its reverse (increasing) ages.
        rates = np.interp(ages, self.age_gy[::-1], self._total_rate()[::-1])
        top = float(rates.max())
        return AgeLaw(tuple(float(a) for a in ages), tuple(round(max(float(r) / top, 0.0), 6) for r in rates))

    def present_over_mean(self):
        """The whole disk's star formation now over its mean over the disk's life (rings weighted by area)."""
        total_rate = self._total_rate()
        return float(total_rate[-1] / np.mean(total_rate))

    def _total_rate(self):
        """The disk's star formation at each step: the rings' rates weighted by their areas (equal widths)."""
        return np.sum(self.radius_scale_lengths[:, None] * self.rate, axis=0)


def infall_time_gy(radius_scale_lengths):
    """
    The infall time scale (Gy) of the rings at `radius_scale_lengths` (array or float): the line of [R5095] in the radius the Milky Way's
    ring of that many scale lengths has, between `tuning.DISK_INFALL_MIN_GY` and `tuning.DISK_INFALL_MAX_GY`.
    """
    radius_kpc = np.asarray(radius_scale_lengths, dtype=float) * tuning.DISK_INFALL_REFERENCE_SCALE_LENGTH_KPC
    line = tuning.DISK_INFALL_SLOPE_GY_PER_KPC * radius_kpc - tuning.DISK_INFALL_OFFSET_GY
    return np.clip(line, tuning.DISK_INFALL_MIN_GY, tuning.DISK_INFALL_MAX_GY)


def disc_speed_kms(stellar_mass_solar):
    """The circular speed (km/s) of a disc of `stellar_mass_solar`: the Milky Way's 220 km/s at its 5e10 solar masses, as the cube root of the mass ([R5096], `M_d` proportional to `V_C^3`)."""
    return tuning.DISK_INFALL_REFERENCE_SPEED_KMS * (stellar_mass_solar / tuning.DISK_INFALL_REFERENCE_MASS_SOLAR) ** (1.0 / 3.0)


def inflow_shapes(radius_scale_lengths, span_gy, steps, speed_kms):
    """
    The gas inflow of each ring over `span_gy` Gy on `steps` equal steps (oldest first), relative to its start: `exp(-t / tau + k t (1 - V / 220))` with
    `infall_time_gy` of the ring, `k` `tuning.DISK_INFALL_MASS_TERM_PER_GY` and `V` the disc's circular speed `speed_kms` ([R5096], equation 21).

    Returns:
        numpy.ndarray: `steps + 1` rows by one column per ring.
    """
    since_start = np.linspace(0.0, span_gy, steps + 1)
    mass_term = tuning.DISK_INFALL_MASS_TERM_PER_GY * (1.0 - speed_kms / tuning.DISK_INFALL_REFERENCE_SPEED_KMS)
    return np.exp(-since_start[:, None] / infall_time_gy(radius_scale_lengths)[None, :] + mass_term * since_start[:, None])


def mass_return_fraction(breaks, slopes):
    """
    The share of a stellar generation's mass that returns to the gas, with instantaneous recycling (Madau and Dickinson
    2014, section 2.1), for a mass function that is a continuous piecewise power law `dN/dm ~ m^-slope` between
    `breaks` (solar masses): stars above `tuning.RECYCLING_MASS_SOL` return all but their remnant, the white dwarf of
    the initial-final relation of `physics.remnant_mass` up to 8 solar masses, a neutron star to `tuning.BLACK_HOLE_RETURN_MASS_SOL`, and a black hole keeps
    its whole mass. (A single Salpeter slope of 2.35 from 0.1 to 100 gives their 0.27.)
    """
    mass = np.geomspace(breaks[0], breaks[-1], 200_001)
    density = imf.density(mass, breaks, slopes)
    remnant = np.where(mass < mass_boundaries.sol("core_collapse_limit"), remnant_mass.white_dwarf_mass_sol(mass),
                       np.where(mass < tuning.BLACK_HOLE_RETURN_MASS_SOL, tuning.NEUTRON_STAR_RETURN_REMNANT_SOL, mass))
    returned = np.where(mass >= tuning.RECYCLING_MASS_SOL, mass - remnant, 0.0)
    log_mass = np.log(mass)
    total = np.trapezoid(mass * mass * density, log_mass)
    return float(np.trapezoid(returned * mass * density, log_mass) / total)


@functools.lru_cache(maxsize=1)
def return_fraction():
    """`mass_return_fraction` for the generator's own mass function (`tuning.IMF_BREAKS_SOL`, `IMF_SLOPES`)."""
    return mass_return_fraction(tuning.IMF_BREAKS_SOL, tuning.IMF_SLOPES)


def _solve_rings(stars_now, inflow, span_gy):
    """
    The rings' histories: for each stellar surface density in `stars_now` (solar masses per pc2), the inflow amplitude
    (the ring's column of `inflow`, rings by steps, times it) for which gas continuity ends with that many stars. Vectorised over rings; the
    amplitude is bracketed by doubling from the least that could do it and bisected.

    Returns:
        tuple: `(rate, gas_now)`, rings by steps and rings.
    """
    steps = inflow.shape[0] - 1
    step_gy = span_gy / steps
    keep = 1.0 - return_fraction()
    coefficient = tuning.STAR_FORMATION_LAW_COEFFICIENT * 1.0e3         # per Gy per pc2, gas in solar masses per pc2
    exponent = tuning.STAR_FORMATION_LAW_EXPONENT
    supply = np.trapezoid(inflow, dx=step_gy, axis=0)
    saturation = tuning.STAR_FORMATION_SATURATION_GAS_MSUN_PC2
    steepening = tuning.STAR_FORMATION_LOW_DENSITY_EXPONENT - exponent

    def thin_gas(gas):
        # The index rises where the gas thins (GEN.450, GEN.999): the law's factor (gas / (gas + saturation)) ^ (m - n).
        return (gas / (gas + saturation)) ** steepening

    def run(amplitude):
        # Semi-implicit: the gas used in a step is taken at the step's end, so a ring so dense that it uses its gas
        # up within a step (a toy galaxy) stays positive and stable.
        gas = np.zeros_like(amplitude)
        rate = np.zeros((steps + 1, len(amplitude)))
        for i in range(steps):
            use = step_gy * keep * coefficient * gas ** (exponent - 1.0) * thin_gas(gas)
            gas = (gas + step_gy * amplitude * inflow[i]) / (1.0 + use)
            rate[i + 1] = coefficient * gas ** exponent * thin_gas(gas)
        return rate, gas

    def formed(rate):
        return keep * np.trapezoid(rate, dx=step_gy, axis=0)

    low = stars_now / supply
    high = low.copy()
    for _ in range(_DOUBLINGS):
        short = formed(run(high)[0]) < stars_now
        if not short.any():
            break
        high = np.where(short, high * 2.0, high)
    low = np.where(high > low, high / 2.0, low)
    for _ in range(_BISECTION_STEPS):
        middle = 0.5 * (low + high)
        under = formed(run(middle)[0]) < stars_now
        low = np.where(under, middle, low)
        high = np.where(under, high, middle)
    rate, gas = run(0.5 * (low + high))
    return rate.T, gas


@functools.lru_cache(maxsize=64)
def history(thin_stellar_mass_solar, disk_scale_length_pc, speed_kms, span_gy=None):
    """
    The thin disk's rings run from `span_gy` ago (`tuning.STAR_FORMATION_AGE_RANGE_GY`'s top by default) to now.

    Args:
        thin_stellar_mass_solar (float): The stars of the thin disk (what the density model integrates to).
        disk_scale_length_pc (float): The stellar disk's exponential scale length.
        speed_kms (float): The disc's circular speed (`disc_speed_kms`), which sets how fast it is fed relative to the Milky Way.
        span_gy (float, optional): How long the disk has been forming stars.

    Returns:
        DiskHistory
    """
    span = float(tuning.STAR_FORMATION_AGE_RANGE_GY[1] if span_gy is None else span_gy)
    rings = tuning.DISK_HISTORY_RINGS
    reach = tuning.DISK_HISTORY_MAX_RADIUS_SCALE_LENGTHS
    x = (np.arange(rings) + 0.5) * reach / rings
    centre = thin_stellar_mass_solar / (2.0 * math.pi * disk_scale_length_pc ** 2)
    stars = centre * np.exp(-x)
    inflow = inflow_shapes(x, span, tuning.DISK_HISTORY_TIME_STEPS, speed_kms)
    rate, gas = _solve_rings(stars, inflow, span)
    return DiskHistory(radius_scale_lengths=x, age_gy=np.linspace(span, 0.0, len(inflow)), rate=rate, gas_now=gas,
                       stars_now=stars)


def stellar_mass_for_scale_length(disk_scale_length_pc):
    """
    The stellar mass of the galaxy whose median disk ([R4048]'s size-mass relation for late types) has this
    scale length, solar masses, held to the mass chain's range. A galaxy with this disk scale length but a size
    drawn off the median differs by a fraction of a dex, which moves the shares by a few percent.
    """
    from planetgen.galaxy import shape_draw
    target = disk_scale_length_pc * shape_draw.EXPONENTIAL_HALF_LIGHT_OVER_SCALE / 1000.0
    low, high = 7.0, 12.5
    for _ in range(60):
        middle = 0.5 * (low + high)
        if shape_draw.size_median_kpc(10.0 ** middle) < target:
            low = middle
        else:
            high = middle
    return 10.0 ** min(max(0.5 * (low + high), 8.2), 11.9)


def stellar_and_halo_masses(shape):
    """`(stellar mass, halo mass)` of the galaxy of this shape, solar masses: its `masses` (from the mass chain) when it
    has them, else those of the median galaxy of its disk scale length. The halo mass is `None` for a galaxy with no
    halo (`mass_chain.has_halo`)."""
    if shape.masses is None:
        stellar_mass = stellar_mass_for_scale_length(shape.disk_scale_length_pc)
        return stellar_mass, mass_chain.median_halo_mass(stellar_mass) if mass_chain.has_halo(stellar_mass) else None
    return shape.masses


def history_of(shape):
    """
    The thin disk's `DiskHistory` for the galaxy of this shape: its `masses` (stellar mass from the mass chain, the halo is not read) when
    it has them, else the median galaxy of its disk scale length; the disk scale length sets the rings, and its bulge and
    thick disk set how much of the stars are the thin disk's. It does not read the density's normalization, so a shape
    and the same shape before `k_norm` was solved share one history (GEN.545 asks for it while `k_norm` is being solved).
    """
    return _history_of(shape._replace(k_norm=1.0))


@functools.lru_cache(maxsize=64)
def _history_of(shape):
    from planetgen.galaxy import density
    length = shape.disk_scale_length_pc
    stellar_mass, _halo_mass = stellar_and_halo_masses(shape)
    masses = density.integrated_masses(shape)
    thin_share = masses["thin_disk"] / (masses["thin_disk"] + masses["thick_disk"] + masses["bulge"])
    return history(stellar_mass * thin_share, length, disc_speed_kms(stellar_mass))


@functools.lru_cache(maxsize=64)
def population_age_laws(shape):
    """
    How each disk population's ages are spread over its window, from the history of the galaxy of this shape (GEN.546):
    `young`, `intermediate` and `old` to their `AgeLaw` (`old` is the thin disc's law mixed with the thick disc's, `thick_ages`, GEN.1025). A population the history does not cover (the bulge, the halo, a
    cluster) has no entry and draws uniformly.

    Args:
        shape (GalaxyShape): As in `population_shares`.

    Returns:
        dict: Population name -> `AgeLaw`.
    """
    from planetgen.galaxy import thick_ages
    found = history_of(shape)
    ages = {**tuning.STELLAR_POPULATION_AGE_RANGES_GY, "old": tuning.THIN_DISC_OLD_WINDOW_GY}
    laws = {name: found.age_law(*ages[name]) for name in tuning.DISK_AGE_LAW_POPULATIONS}
    laws["old"] = thick_ages.old_age_law(shape, laws["old"])    # GEN.1025: the thick disc's stars reach 13 Gy
    return laws


@functools.lru_cache(maxsize=64)
def population_shares(shape):
    """
    The share of the thin disk's stars in each disk population (`young`, `intermediate`, `old`), from the history of the
    galaxy of this shape (`history_of`). They add up to 1.

    Args:
        shape (GalaxyShape): The galaxy's shape.

    Returns:
        dict: Population name -> share.
    """
    found = history_of(shape)
    ages = {**tuning.STELLAR_POPULATION_AGE_RANGES_GY, "old": tuning.THIN_DISC_OLD_WINDOW_GY}    # GEN.1025: the thin disc's own window
    shares = {name: found.share_between(*ages[name]) for name in tuning.STELLAR_POPULATION_ARM_AMPLITUDE}
    total = sum(shares.values())
    return {name: value / total for name, value in shares.items()}
