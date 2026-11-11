# planetgen/physics/disk_size.py

"""
The outer edge of a star's protoplanetary gas disk (GEN.1520).

The edge is the cut-off radius of the gas disk the planets, belts and comets of a star formed in. It follows the star's mass at formation, not its snow line:

    R_edge = `EDGE_AT_SUN_AU` x (M / Msun) ^ `EDGE_MASS_EXPONENT`

and each star draws its own edge around that median, log-normal, `EDGE_SCATTER_DEX` dex wide (a factor of 2 either way) and cut at `EDGE_SCATTER_MAX_SIGMA` sigma. The draw is
made once per star from a stream keyed by the star (no stored column: any caller gets the same edge), so it moves no other draw.

Evidence (`reports/gen1520-disk-size-by-mass-2026-11-11.md`): the dust radius of ALMA disks is 27 to 60 au at the median (Hendler et al. 2020), the CO gas disk is a median 2.5 times
the dust disk (Sanchis et al. 2021), and the size rises with star mass as M^0.2 to 0.6 (Hendler 2020 per region, Andrews et al. 2018: 0.6); the planet-formation models of Burn et al.
(2021) cut the gas disk at 62 au for a Sun-like star and 15 au at 0.1 Msun (M^0.625). 100 au is the geometric middle of the 62 au model and the 200 au or so the observed
gas radii give, and the exponent 0.5 the middle of the direct fits.
"""

from planetgen.physics import mass_loss, unit_registry as units
from planetgen.util import draw

_kg_to_msun = units.converter("kg", "Msun")

EDGE_AT_SUN_AU = 100.0
"""float: The median gas-disk edge of a 1 solar-mass star, AU. [Stand-in DS1, tier Marginal Unverified, tweak L1: between the 62 au of Burn et al. 2021 and the roughly 200 au the observed CO radii of 1-3 Myr disks give.]"""

EDGE_MASS_EXPONENT = 0.5
"""float: The edge rises as the star's mass to this power. [Stand-in DS2, tier Confident Unverified, tweak L1: the direct fits of dust radius against mass are 0.2 to 0.57 (Hendler 2020), 0.6 (Andrews 2018); the models use 0.625.]"""

EDGE_SCATTER_DEX = 0.3
"""float: The scatter of a star's own edge around the median, dex in log10: a factor of 2 either way. [Stand-in DS3, tier Marginal Unverified, tweak L1: real disks of one mass vary by a factor of about 4 either way in radius (Hendler 2020's scatter 0.14 to 0.43 dex in the fits).]"""

EDGE_SCATTER_MAX_SIGMA = 2.0
"""float: The draw is cut at this many sigma (a factor of 4 either way). [Stand-in DS4, tier Marginal Unverified, tweak L0.]"""


def median_edge_au(mass_solar):
    """The median gas-disk edge, AU, of a star of `mass_solar` solar masses at formation."""
    if not mass_solar > 0.0:
        raise ValueError("the star's mass must be positive")
    return EDGE_AT_SUN_AU * mass_solar ** EDGE_MASS_EXPONENT


def edge_au_from(mass_solar, z):
    """The edge, AU, of a star of `mass_solar` that drew `z` sigma from the median."""
    return median_edge_au(mass_solar) * 10.0 ** (EDGE_SCATTER_DEX * z)


def draw_edge_au(mass_solar, stream):
    """One star's edge, AU: the median of `median_edge_au` times 10 to a normal draw of `EDGE_SCATTER_DEX` dex from `stream`, redrawn while it lies beyond `EDGE_SCATTER_MAX_SIGMA` sigma."""
    while True:
        z = stream.gauss(0.0, 1.0)
        if abs(z) <= EDGE_SCATTER_MAX_SIGMA:
            return edge_au_from(mass_solar, z)
        draw.redrawn("disk_edge.scatter")


def _members(star):
    return getattr(star, "stars", None) or [star]


def formation_mass_solar(star):
    """The mass, solar masses, `star` (a star or a close pair's proxy) had when its disk formed: its members' mass at formation, summed (a white dwarf is counted at its progenitor's)."""
    return _kg_to_msun(sum(mass_loss.member_formation_mass_kg(member) for member in _members(star)))


def edge_au(star):
    """
    The outer edge, AU, of the gas disk of `star` (a star, a compact remnant or a close pair's proxy): `draw_edge_au` for its mass at formation, from a stream keyed by its
    members' mass, age and place in the galaxy, so every caller gets the same edge and no run's draw moves (the way the star's own snow line is drawn).
    """
    key = "|".join(f"{member.mass!r}:{member.age!r}:{member.galactic_center_dist_ly!r}" for member in _members(star))
    return draw_edge_au(formation_mass_solar(star), draw.Stream(f"disk_edge:{key}"))
