# planetgen/galaxy/scatter_limits.py

"""
The scatter limits by place (DB.54, docs/design/realistic-galactic-core.md 5.1)
==============================================================================

The galaxy's star scatter draws its line at one luminosity floor and one mass limit (the baseline,
`galaxy_shape.bright_star_min_luminosity_sol` and `bright_star_mass_limit_sol`). Three more settings of the galaxy move
that line by place: a luminosity floor and a mass limit for the core and a luminosity floor for the halo (the halo has
no mass limit of its own). The core's two are `None` until set, and then follow the baseline times
`tuning.CORE_TIER_FACTOR`; the halo's is a number (default `tuning.HALO_MIN_LUMINOSITY_SOL`, never under
`tuning.HALO_MIN_LUMINOSITY_RANGE_SOL[0]`).

Which place takes which limit (GEN.333)
=======================================

The scatter walks the galaxy a *band* at a time (a shell and a band index: the unit the object-first sampler bounds
its density over), so a band is the unit of a tier. `tier_indices` gives each band of a shell its tier from the
population mix at the band's centre line, averaged over azimuth: the **core** tier when the two nuclear terms are at
least `tuning.TIER_DOMINANCE_SHARE` of the band's density, the **halo** tier when the halo is, the **baseline** tier
everywhere else (the thin and thick disc, the bulge, the arms). One tier per band is a pure function of the address, so
the scatter's passes, the backfill, the fill and the sector page all read the same limits from it with no state kept.

`ScatterLimits` holds the baseline the scatter ran at (its luminosity floor and mass limit, `None` before any scatter
ran) and the three settings of `galaxy_shape`, and gives the limits of a tier and of an address. The core takes its
own floor and mass limit, the halo its own floor and the baseline's mass limit.
"""

import functools
import math
from dataclasses import dataclass

import numpy as np

from planetgen import tuning
from planetgen.galaxy import density as galaxyDensity
from planetgen.galaxy.component_field import POPULATIONS_IN_THE_GATE
from planetgen.galaxy.spherical import band_center_pc, wedge_count

CORE, BASELINE, HALO = 0, 1, 2
"""int: The tier indices, in the order the three settings are listed."""

TIER_NAMES = ("core", "baseline", "halo")

TIER_AZIMUTH_BINS = 8
"""int: Longitude bins a band's population mix is averaged over to choose its tier (fewer where the band has
fewer wedges). The tier only needs the dominant term, so it is coarser than the scatter's own density bins."""


def baseline_min_luminosity_sol(level_sol):
    """The baseline luminosity floor: the scatter's level when one has run, else the tuning's default."""
    return tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL if level_sol is None else float(level_sol)


def baseline_mass_limit_sol(limit_sol):
    """The baseline mass limit: the scatter's when it had a mass pass, else the tuning's default."""
    return tuning.PHENOMENON_MIN_MASS_SOLAR if limit_sol is None else float(limit_sol)


def core_min_luminosity_sol(skeleton, level_sol):
    """The core's luminosity floor: the one set, else the baseline (`level_sol`, the scatter's level or `None`) times
    `tuning.CORE_TIER_FACTOR`."""
    if skeleton.core_min_luminosity_sol is not None:
        return float(skeleton.core_min_luminosity_sol)
    return baseline_min_luminosity_sol(level_sol) * tuning.CORE_TIER_FACTOR


def core_mass_limit_sol(skeleton, limit_sol):
    """The core's mass limit: the one set, else the baseline (`limit_sol`, the scatter's or `None`) times
    `tuning.CORE_TIER_FACTOR`."""
    if skeleton.core_mass_limit_sol is not None:
        return float(skeleton.core_mass_limit_sol)
    return baseline_mass_limit_sol(limit_sol) * tuning.CORE_TIER_FACTOR


def halo_min_luminosity_sol(skeleton):
    """The halo's luminosity floor."""
    return float(skeleton.halo_min_luminosity_sol)


def check_halo_min_luminosity(value):
    """Raises `ValueError` naming the reason when `value` is not a number inside the halo floor's range."""
    low, high = tuning.HALO_MIN_LUMINOSITY_RANGE_SOL
    try:
        ok = low <= float(value) <= high
    except (TypeError, ValueError):
        ok = False
    if not ok:
        raise ValueError(f"the halo luminosity floor must be from {low:g} to {high:g} solar luminosities, "
                         f"got {value!r}")


_NUCLEAR = ("nuclear_cluster", "nuclear_disc")


@functools.lru_cache(maxsize=TIER_AZIMUTH_BINS + 1)
def _azimuths(count):
    """`(cos, sin)` of the `count` longitudes `(k + 0.5)` of a turn a band's mix is read at, as arrays; the very floats
    the plain loop over `math.cos` and `math.sin` gave (numpy's own trigonometry can differ in the last bit)."""
    thetas = [(k + 0.5) * 2.0 * math.pi / count for k in range(count)]
    return (np.array([math.cos(theta) for theta in thetas]), np.array([math.sin(theta) for theta in thetas]))


def _tiers_of_bands(shape, shell_index, bands, edge_pc):
    """The tier index of each of `bands` (a list of band indices) of a shell, as an array."""
    return _tiers_of_places(shape, [(shell_index, band) for band in bands], edge_pc)


def _tiers_of_places(shape, places, edge_pc):
    """The tier index of each of `places` (a list of `(shell, band)`), as an array: one pass over the density
    arrays however many shells the places lie in."""
    if not places:
        return np.zeros(0, dtype=np.int8)
    counts = np.array([min(wedge_count(shell_index, band), TIER_AZIMUTH_BINS) for shell_index, band in places])
    centres = [band_center_pc(shell_index, band, edge_pc) for shell_index, band in places]
    planes = np.repeat(np.array([plane for plane, _z in centres]), counts)
    zs = np.repeat(np.array([z for _plane, z in centres]), counts)
    cosines = np.concatenate([_azimuths(int(count))[0] for count in counts])
    sines = np.concatenate([_azimuths(int(count))[1] for count in counts])
    xs, ys = planes * cosines, planes * sines
    found = galaxyDensity.population_densities_array(xs, ys, zs, shape)
    kept = {name: np.maximum(array, 0.0) for name, array in found.items()}
    for kind in POPULATIONS_IN_THE_GATE:
        kept[kind] = np.zeros_like(kept[kind])
    offsets = np.concatenate(([0], np.cumsum(counts)[:-1]))

    def per_band(array):
        return np.add.reduceat(array, offsets)

    return _dominant_tiers({name: per_band(array) for name, array in kept.items()})


def _dominant_tiers(summed):
    """The tier index of each place from its populations' densities (`{population: array}`, a gate kind's counted
    as none and a negative term as none): the nuclear terms' share, then the halo's, against `TIER_DOMINANCE_SHARE`."""
    total = sum(array for name, array in summed.items() if name not in POPULATIONS_IN_THE_GATE)
    nuclear = sum(summed[name] for name in _NUCLEAR)
    halo = summed["halo"]
    wanted = tuning.TIER_DOMINANCE_SHARE * total
    tiers = np.full(np.shape(total), BASELINE, dtype=np.int8)
    tiers[(total > 0.0) & (halo >= wanted)] = HALO
    tiers[(total > 0.0) & (nuclear >= wanted)] = CORE
    return tiers


REACH_AZIMUTH_BINS = 2
"""int: Longitudes a group of bands beyond the outline is probed at (`probe_densities`): far from the disk the
density is nearly round, and a group is bounded, not measured."""


def probe_densities(shape, plane, z, azimuths=REACH_AZIMUTH_BINS):
    """
    The population densities at each of the places `(plane, z)` (arrays of cylindrical radius and height, parsecs),
    each averaged over `azimuths` evenly spaced longitudes, and the tier those densities give the place
    (`tier_indices`' rule): `(means, tiers)`, `means` `{population: array}` with a gate kind's and a negative term
    counted as none. What a group of bands beyond the outline is bounded and weighed from (GEN.322).
    """
    plane, z = np.asarray(plane, dtype=float), np.asarray(z, dtype=float)
    theta = (np.arange(azimuths) + 0.5) * 2.0 * math.pi / azimuths
    found = galaxyDensity.population_densities_array(
        (plane[:, None] * np.cos(theta)).ravel(), (plane[:, None] * np.sin(theta)).ravel(),
        np.repeat(z, azimuths), shape)
    means = {}
    for name, array in found.items():
        array = np.zeros_like(array) if name in POPULATIONS_IN_THE_GATE else np.maximum(array, 0.0)
        means[name] = array.reshape(len(plane), azimuths).mean(axis=1)
    return means, _dominant_tiers(means)


@functools.lru_cache(maxsize=4096)
def _shell_tiers(shape, edge_pc, shell_index, high_band):
    tiers = _tiers_of_bands(shape, shell_index, list(range(-high_band, high_band + 1)), edge_pc)
    tiers.flags.writeable = False
    return tiers


@functools.lru_cache(maxsize=65536)
def _band_tier(shape, edge_pc, shell_index, band_index):
    return int(_tiers_of_bands(shape, shell_index, [band_index], edge_pc)[0])


def tier_indices(shape, shell_index, high_band, edge_pc):
    """The tier index (`CORE`, `BASELINE` or `HALO`) of each band `-high_band .. high_band` of a shell (index
    `band + high_band`), an array."""
    return _shell_tiers(shape, edge_pc, shell_index, high_band)


def tiers_of(shape, shell_index, bands, edge_pc):
    """The tier index of each of `bands` of a shell (a list of band indices), an array: what `tier_indices` gives those
    bands, without working out the others."""
    return _tiers_of_bands(shape, shell_index, list(bands), edge_pc)


def tiers_of_places(shape, places, edge_pc):
    """The tier index of each of `places` (a list of `(shell, band)`, shells of any number), an array: what `tier_index`
    gives each, worked out together."""
    return _tiers_of_places(shape, list(places), edge_pc)


def tier_index(shape, shell_index, band_index, edge_pc):
    """The tier index of one band."""
    return _band_tier(shape, edge_pc, shell_index, band_index)


@dataclass(frozen=True)
class TierLimit:
    """
    What one tier's scatter draws its line at.

    Attributes:
        min_luminosity_sol (float): Every star at least this bright is pre-placed.
        mass_limit_sol (float or None): Every star born at this mass or more is pre-placed, whatever its luminosity;
            `None` when the scatter had no mass pass.
    """
    min_luminosity_sol: float
    mass_limit_sol: float | None


@dataclass(frozen=True)
class ScatterLimits:
    """
    The scatter's limits by place for one galaxy: the baseline it ran at and the three tier limits.

    Attributes:
        shape (GalaxyShape): The galaxy's shape, which the tiers are chosen from.
        edge_pc (float): The sector edge.
        tiers (tuple): `(core, baseline, halo)` `TierLimit`s, or `None` before any scatter ran (nothing is
            pre-placed, so no limit applies).
    """
    shape: object
    edge_pc: float
    tiers: tuple | None

    @classmethod
    def of(cls, skeleton, level_sol, mass_limit_sol):
        """
        The limits for `skeleton` when the scatter ran (or is to run) at baseline floor `level_sol` and mass limit
        `mass_limit_sol` (`None` for no mass pass). `level_sol` `None` means no scatter ran: no tier has a limit.
        """
        if level_sol is None:
            return cls(skeleton.shape, skeleton.edge_pc, None)
        level = float(level_sol)
        if mass_limit_sol is None:
            core_mass = halo_mass = baseline_mass = None
        else:
            baseline_mass = float(mass_limit_sol)
            core_mass = core_mass_limit_sol(skeleton, baseline_mass)
            halo_mass = baseline_mass
        return cls(skeleton.shape, skeleton.edge_pc, (
            TierLimit(core_min_luminosity_sol(skeleton, level), core_mass),
            TierLimit(level, baseline_mass),
            TierLimit(halo_min_luminosity_sol(skeleton), halo_mass)))

    def limit(self, tier):
        """The `TierLimit` of a tier index, `None` before any scatter ran."""
        return None if self.tiers is None else self.tiers[tier]

    def tiers_in(self, shell_index, high_band):
        """`tier_indices` for a shell of this galaxy."""
        return tier_indices(self.shape, shell_index, high_band, self.edge_pc)

    def tiers_of(self, shell_index, bands):
        """`tiers_of` for a shell of this galaxy."""
        return tiers_of(self.shape, shell_index, bands, self.edge_pc)

    def tiers_of_places(self, places):
        """`tiers_of_places` for this galaxy."""
        return tiers_of_places(self.shape, places, self.edge_pc)

    def tier_at(self, shell_index, band_index):
        """The tier index of one band."""
        return tier_index(self.shape, shell_index, band_index, self.edge_pc)

    def limit_at(self, shell_index, band_index):
        """The `TierLimit` of one band, `None` before any scatter ran. A place with no band (`None`) takes the
        baseline's."""
        if self.tiers is None:
            return None
        if shell_index is None or band_index is None:
            return self.tiers[BASELINE]
        return self.limit(self.tier_at(shell_index, band_index))

    def level_at(self, shell_index, band_index):
        """The luminosity floor the scatter drew a band's stars down to, `None` before any scatter ran."""
        found = self.limit_at(shell_index, band_index)
        return None if found is None else found.min_luminosity_sol

    def mass_limit_at(self, shell_index, band_index):
        """The mass limit of a band's mass pass, `None` when it had none (or no scatter ran)."""
        found = self.limit_at(shell_index, band_index)
        return None if found is None else found.mass_limit_sol


def tier_of_densities(densities):
    """
    The tier index of a band from the population densities at its centre line (a list of `{population: density}`,
    one per longitude bin, a gate kind's counted as none and a negative term as none, as the scatter reads them): the
    same rule `tier_indices` applies to a whole shell at once.
    """
    total = nuclear = halo = 0.0
    for found in densities:
        for name, value in found.items():
            if name in POPULATIONS_IN_THE_GATE:
                continue
            value = max(value, 0.0)
            total += value
            if name in _NUCLEAR:
                nuclear += value
            elif name == "halo":
                halo += value
    if not total > 0.0:
        return BASELINE
    wanted = tuning.TIER_DOMINANCE_SHARE * total
    if nuclear >= wanted:
        return CORE
    if halo >= wanted:
        return HALO
    return BASELINE


@dataclass(frozen=True)
class PassBand:
    """
    What one tier's part of a scatter pass draws: the stars at least `min_luminosity_sol` bright and (when
    `max_luminosity_sol` is not `None`) dimmer than it, born in `mass_range` (`(low, high)` solar masses, an end
    `None` for no limit; `None` for any).
    """
    min_luminosity_sol: float
    max_luminosity_sol: float | None
    mass_range: tuple | None


@dataclass(frozen=True)
class TieredPass:
    """
    A scatter pass by place: `limits` chooses each band's tier and `bands` holds the `PassBand` of the core,
    baseline and halo tiers, a tier with `None` drawing nothing in this pass.
    """
    limits: ScatterLimits
    bands: tuple

    def tiers_in(self, shell_index, high_band):
        return self.limits.tiers_in(shell_index, high_band)


def mass_pass(limits, min_luminosity_sol):
    """The mass pass: every star born at or above its tier's mass limit, whatever its luminosity (at least
    `min_luminosity_sol`, the dimmest the star models make). A tier with no mass limit draws nothing."""
    return TieredPass(limits, tuple(
        None if limits.tiers[tier].mass_limit_sol is None
        else PassBand(min_luminosity_sol, None, (limits.tiers[tier].mass_limit_sol, None))
        for tier in (CORE, BASELINE, HALO)))


def luminosity_pass(limits):
    """The luminosity pass: every star at or above its tier's floor that was born under its tier's mass limit."""
    return TieredPass(limits, tuple(
        PassBand(found.min_luminosity_sol, None,
                 None if found.mass_limit_sol is None else (None, found.mass_limit_sol))
        for found in limits.tiers))


def band_pass(newer, older):
    """
    The band a staged scatter adds (`planetgen scatter --band --down-to`): per tier, the stars from its new floor
    up to, not including, its old one. A tier whose floor did not fall (the halo's, which does not follow the
    baseline; a core floor that was set) draws nothing.
    """
    bands = []
    for new, old in zip(newer.tiers, older.tiers):
        if new.min_luminosity_sol >= old.min_luminosity_sol:
            bands.append(None)
        else:
            bands.append(PassBand(new.min_luminosity_sol, old.min_luminosity_sol,
                                  None if new.mass_limit_sol is None else (None, new.mass_limit_sol)))
    return TieredPass(newer, tuple(bands))


def tiered_to_json(tiered):
    """A `TieredPass` as plain JSON (a unit on the work queue, PERF.97): the tier limits and the three bands; the
    galaxy's shape is not in it, the unit reads that from the database (`tiered_from_json`). `None` for `None`."""
    if tiered is None:
        return None
    tiers = None if tiered.limits.tiers is None else [[t.min_luminosity_sol, t.mass_limit_sol]
                                                       for t in tiered.limits.tiers]
    bands = [None if band is None else
             [band.min_luminosity_sol, band.max_luminosity_sol,
              None if band.mass_range is None else list(band.mass_range)] for band in tiered.bands]
    return {"tiers": tiers, "bands": bands}


def tiered_from_json(data, skeleton):
    """The `TieredPass` `tiered_to_json` wrote, for the galaxy `skeleton` (its shape and sector edge)."""
    if data is None:
        return None
    tiers = None if data["tiers"] is None else tuple(TierLimit(low, mass) for low, mass in data["tiers"])
    bands = tuple(None if band is None else
                  PassBand(band[0], band[1], None if band[2] is None else tuple(band[2])) for band in data["bands"])
    return TieredPass(ScatterLimits(skeleton.shape, skeleton.edge_pc, tiers), bands)
