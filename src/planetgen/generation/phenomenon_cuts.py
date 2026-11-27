# planetgen/generation/phenomenon_cuts.py

"""
The phenomenon scatter's mass cut by place (GEN.377).

The scatter places the neutron stars and black holes at or above a mass cut and leaves the share below it to the sector's fill
(`phenomenon_scatter.below_cut_draws`). It kept one cut for the whole galaxy while GEN.333 gave the stars three limits by place:
the star scatter draws the core, the baseline and the halo each at its own limit, so a core band with a heavier limit had its stars
placed further up the mass scale than a disk band did, and the compact objects' cut stayed the same everywhere. `Cuts` is the same
three limits for the compact objects, from the same settings the Generate page shows (the baseline cut, `--core-mass-limit` and the
halo, which has no mass limit of its own and takes the baseline's: `scatter_limits.ScatterLimits.of`), with each band's tier chosen
by the dominant population of the band as for stars (`scatter_limits.tier_indices`), so one band takes one cut.

A galaxy whose three cuts are equal (`uniform`) draws exactly what one cut drew before, with no tier looked up.

The walk by band (`phenomenon_scatter._band_rows`) gives each band its tier's rate. The two scatters that place candidates and keep
them (the reach beyond the outline, `phenomenon_reach`, and the clusters and satellites, `component_pass`) draw the candidates at the
lowest cut's rate and keep each with the probability its tier's share over that (`keep_share`), as the star pass keeps by tier.
"""

from planetgen import tuning
from planetgen.galaxy import scatter_limits


class Cuts:
    """
    The lowest mass of a neutron star or black hole the scatter places, by tier: `tiers` is `(core, baseline, halo)` in solar
    masses, and `limits` the galaxy's `scatter_limits.ScatterLimits` that gives a band its tier (`None` for one cut everywhere).
    """

    __slots__ = ("tiers", "limits")

    def __init__(self, tiers, limits=None):
        self.tiers = tuple(float(value) for value in tiers)
        self.limits = limits

    @classmethod
    def uniform(cls, cut_solar):
        """One cut for the whole galaxy."""
        cut = float(cut_solar)
        return cls((cut, cut, cut))

    @classmethod
    def of(cls, skeleton, cut_solar):
        """
        The cuts of `skeleton`'s galaxy when the scatter's baseline cut is `cut_solar`: the baseline's, the halo's (the same) and
        the core's (the setting the galaxy has, else the baseline times `tuning.CORE_TIER_FACTOR`: `ScatterLimits.of`).
        """
        found = scatter_limits.ScatterLimits.of(skeleton, tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL, float(cut_solar))
        return cls([tier.mass_limit_sol for tier in found.tiers], found)

    @classmethod
    def from_tiers(cls, tiers):
        """The cuts of a record that holds only the three numbers (a log line, a plan's detail): they name no band's tier."""
        return cls(tiers)

    @property
    def is_uniform(self):
        """Whether every tier has the same cut: one cut is then used for every band and no tier is looked up."""
        return len(set(self.tiers)) == 1

    @property
    def baseline(self):
        """The baseline tier's cut: the number the scatter's settings record."""
        return self.tiers[scatter_limits.BASELINE]

    @property
    def highest(self):
        """The highest of the three cuts: the one that leaves the most below it, which the promoted compact pairs' candidate counts are drawn at (GEN.263)."""
        return max(self.tiers)

    @property
    def lowest(self):
        """The lowest of the three cuts: the one that places the most, which the candidate counts and the bounds are drawn at."""
        return min(self.tiers)

    def tiers_in(self, shell_index, high_band):
        """The tier index of each band `-high_band .. high_band` of a shell (index `band + high_band`), an array."""
        return self.limits.tiers_in(shell_index, high_band)

    def tiers_of_places(self, places):
        """The tier index of each of `places` (a list of `(shell, band)`), an array."""
        return self.limits.tiers_of_places(places)

    def cut_at(self, shell_index, band_index):
        """The cut of one band; a place with no band takes the baseline's."""
        if self.is_uniform:
            return self.tiers[0]
        if shell_index is None or band_index is None:
            return self.baseline
        return self.tiers[self.limits.tier_at(shell_index, band_index)]

    def cut_of_tier(self, tier):
        """The cut of a tier index."""
        return self.tiers[tier]


def as_cuts(value):
    """`value` as `Cuts`: a number is one cut for the whole galaxy, `None` the tuning's default, a `Cuts` itself."""
    if isinstance(value, Cuts):
        return value
    return Cuts.uniform(tuning.PHENOMENON_MIN_MASS_SOLAR if value is None else value)
