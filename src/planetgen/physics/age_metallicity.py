# planetgen/physics/age_metallicity.py

"""
The age-metallicity relation of the disc, by branch (GEN.744 part 2, GEN.944 part 1).

A star's [Fe/H] was drawn from its population alone (`metallicity.sample_fe_h`), independent of its age inside the population's window.
In the solar neighbourhood the two are tied, and differently on the two branches (Haywood et al. 2013, A&A 560, A109, Section 3.3 and
Fig. 9):

- **The thin disc** (`stellar_abundances.BRANCH_THIN`) gains metals slowly, 0.025 dex per Gyr, with a wide scatter.
- **The thick disc** (`BRANCH_HIGH_ALPHA`, older than `tuning.STELLAR_THICK_START_GY`, 8 Gyr) shows a tight relation that rises 0.15
  dex per Gyr over the 4 to 5 Gyr it took to form, "also implying a decrease by a factor of 5-6 in the production of iron after 8 Gyr".

The relation is the mean `mean_fe_h(age, branch)`: a line through the branch's anchor, at the slope the paper gives. The anchors (the thin
disc's at birth, the thick disc's at 8 Gyr) and the scatters are stand-ins (AR1 to AR3) chosen to match the populations the generator has
today (`tuning.METALLICITY_BY_POPULATION`: 0 for the young and intermediate disc, -0.2 for the old) and the thin disc's lower limit near -0.7
(Bensby et al. 2014). `sample_fe_h` draws about it. Nothing in generation calls this yet: the wiring draws the age first, then the [Fe/H]
(the order of GEN.553 reverses for the disc), which changes every disc star's stream and so rides train 3.
"""

from planetgen import tuning
from planetgen.physics import stellar_abundances

THIN_SLOPE_DEX_PER_GY = -0.025
"""float: d[Fe/H]/d(age) of the thin disc, dex per Gyr (Haywood et al. 2013: 0.025 dex per Gyr of metal gain, so older is poorer)."""
THICK_SLOPE_DEX_PER_GY = -0.15
"""float: d[Fe/H]/d(age) of the thick disc (Haywood et al. 2013: about 0.15 dex per Gyr)."""
THIN_BIRTH_FE_H = 0.05
"""float: The thin disc's mean [Fe/H] at age 0. [Stand-in AR1, tier T6, tweak L1: the young and intermediate disc of today average 0.]"""
THICK_ANCHOR_FE_H = -0.10
"""float: The thick disc's mean [Fe/H] at `tuning.STELLAR_THICK_START_GY`, where it ends. [Stand-in AR2, tier T6, tweak L1: Haywood's Fig. 9
runs from about -0.8 at 12 Gyr to about 0 at 8 Gyr.]"""
SIGMA = {stellar_abundances.BRANCH_THIN: 0.15, stellar_abundances.BRANCH_HIGH_ALPHA: 0.10}
"""dict: The scatter about the mean, dex, by branch. [Stand-in AR3, tier T6, tweak L1: the thin disc's is the generator's own 0.15; the
thick disc's is tighter, "a much tighter correlation".]"""
THICK_OLDEST_GY = 13.0
"""float: The oldest thick-disc age; the relation is held at its end beyond it."""


def mean_fe_h(age_gy, branch):
    """The mean [Fe/H] of the disc at this age on this branch."""
    if branch == stellar_abundances.BRANCH_THIN:
        return THIN_BIRTH_FE_H + THIN_SLOPE_DEX_PER_GY * age_gy
    if branch == stellar_abundances.BRANCH_HIGH_ALPHA:
        age = min(max(age_gy, tuning.STELLAR_THICK_START_GY), THICK_OLDEST_GY)
        return THICK_ANCHOR_FE_H + THICK_SLOPE_DEX_PER_GY * (age - tuning.STELLAR_THICK_START_GY)
    raise ValueError(f"unknown branch {branch!r}; expected one of {stellar_abundances.BRANCHES}")


def sample_fe_h(age_gy, branch, normal):
    """A disc star's [Fe/H]: the mean at its age plus `normal` (a standard normal draw) times the branch's scatter, held in
    `tuning.METALLICITY_RANGE_FE_H`."""
    low, high = tuning.METALLICITY_RANGE_FE_H
    return min(max(mean_fe_h(age_gy, branch) + SIGMA[branch] * normal, low), high)
