# planetgen/physics/roche.py

"""
Roche geometry (GEN.269): the [R1222] Roche-lobe radius and the
filling factor of a star in a pair. The interaction classifier of
`planetgen.generation.multistar` reads them: a star whose radius reaches its
lobe at the pair's closest approach feeds its partner (semi-detached, or
contact when both reach theirs).
"""

import math


def roche_lobe_fraction(mass_self, mass_other):
    """The Roche-lobe radius of the first body as a fraction of the separation [R1222]:
    `0.49 q^(2/3) / (0.6 q^(2/3) + ln(1 + q^(1/3)))` with `q = mass_self / mass_other`."""
    q13 = (mass_self / mass_other) ** (1.0 / 3.0)
    q23 = q13 * q13
    return 0.49 * q23 / (0.6 * q23 + math.log(1.0 + q13))


def roche_lobe_radius_au(separation_au, mass_self, mass_other):
    """The Roche-lobe radius in AU of the first body at this separation."""
    return separation_au * roche_lobe_fraction(mass_self, mass_other)


def filling_factor(radius_au, separation_au, mass_self, mass_other):
    """The star's radius over its Roche-lobe radius: 1 or more means it fills its lobe."""
    return radius_au / roche_lobe_radius_au(separation_au, mass_self, mass_other)


def separation_for_filling(radius_au, filling, mass_self, mass_other):
    """The separation at which the star has this filling factor."""
    return radius_au / (filling * roche_lobe_fraction(mass_self, mass_other))
