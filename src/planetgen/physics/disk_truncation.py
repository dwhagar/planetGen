# planetgen/physics/disk_truncation.py

"""
How far out a companion star lets a disk reach (GEN.1451)
=========================================================

A star with another star on a bound orbit holds a disk (and so the belts that grew in it) only inside the radius where the companion's tide takes over. The planets' stability
window (`planet_zones`, Holman and Wiegert, stored as `Star.a_crit_au`) says how far a body can *stay*; the disk the belts are made of ends nearer, where the gas and the
solids it carries are cut off by the companion at its closest approach. A belt's stretch is cut there before its mass is worked out (`belt_origin.build`), so a belt inside a close
companion holds the solids of the cut stretch, not of the whole disk.

**The law.** The truncation radius is a share `TRUNCATION_OF_ROCHE_LOBE` of the Roche lobe of the host at the companion's periastron, with the lobe from the fit of [R1222],
`r_L / a = 0.49 q^(2/3) / (0.6 q^(2/3) + ln(1 + q^(1/3)))`, `q` the host's mass over the companion's and `a` the periastron separation `axis (1 - e)`. A disk round one star of an equal
pair then ends near `0.28` of the pair's periastron, the third or so that disks of young binaries are seen cut at. The tightest companion sets the radius; a lone star, or a star the
system has attached no companion to (`radiation_field.companions_of`), has no truncation (`math.inf`).

What it leaves out (stated; stand-in rows in the gap plan): the companion's orbit is taken as it stands today (not as it was while the disk lasted); a disk and its companion are taken
coplanar (a mutual inclination only moves the cut further out); the factor does not depend on the disk's viscosity or Reynolds number; a subtree of stars is one point (as
`radiation_field` takes it).
"""

import math

from planetgen.physics import radiation_field

TRUNCATION_OF_ROCHE_LOBE = 0.733
"""float: The share of the host's Roche lobe at periastron where a circumstellar disk ends (disks of binaries truncate at about 0.7 of the Roche lobe in simulations of the disk's tidal
cut-off, read second-hand). [Stand-in DT1, tier Marginal unverified, tweak L0.]"""


def roche_lobe_share(mass_ratio):
    """The Roche lobe of a star over the separation, for `mass_ratio` its mass over its companion's [R1222]."""
    if not mass_ratio > 0.0:
        raise ValueError("the mass ratio must be positive")
    third = mass_ratio ** (1.0 / 3.0)
    return 0.49 * third * third / (0.6 * third * third + math.log(1.0 + third))


def truncation_radius_au(host_mass_kg, companion_mass_kg, axis_au, eccentricity):
    """The radius, AU, where a companion of `companion_mass_kg` on an orbit of semi-major axis `axis_au` and `eccentricity` round the host cuts the host's disk."""
    if not (host_mass_kg > 0.0 and companion_mass_kg > 0.0 and axis_au > 0.0 and 0.0 <= eccentricity < 1.0):
        return math.inf
    periastron = axis_au * (1.0 - eccentricity)
    return TRUNCATION_OF_ROCHE_LOBE * roche_lobe_share(host_mass_kg / companion_mass_kg) * periastron


def star_truncation_au(star):
    """The radius where the companions that `radiation_field` attached to `star` cut its disk, AU: the tightest of them, `math.inf` for a star with none."""
    radius = math.inf
    for companion in radiation_field.companions_of(star):
        if companion.dark:      # a glowing dark member of the tree (GEN.1501) lights the planets; the gravity of every member is the work of GEN.1446 to GEN.1462
            continue
        mass = getattr(companion.star, "mass", None)
        host = getattr(star, "mass", None)
        if mass and host:
            radius = min(radius, truncation_radius_au(host, mass, companion.axis_au, companion.eccentricity))
    return radius


def star_disk_limit_au(star):
    """Where the companions of `star` end a belt's disk, AU: the radius they cut the disk at (`star_truncation_au`) or the limit past which no orbit holds (`Star.a_crit_au`, Holman and
    Wiegert), whichever is nearer, so the mass of a belt is that of the stretch that is left in either case. `math.inf` for a lone star (it has no `a_crit_au`)."""
    stable = getattr(star, "a_crit_au", None)
    return min(star_truncation_au(star), math.inf if stable is None else stable)
