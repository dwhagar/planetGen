# planetgen/physics/white_dwarf_merger.py

"""
The double white dwarf merger products among the white dwarfs (GEN.708)
========================================================================

A white dwarf of 0.9 solar masses or more is as likely a merger of two white dwarfs as a single star's remnant, and the share
depends on the mass: [R5355] finds 25 percent of the white dwarfs of 0.9 to 1 solar mass to be merger products, 36 percent at 1
to 1.1, 78 percent at 1.1 to 1.2, 49 percent at 1.2 to 1.3 and 38 percent above (`tuning.WHITE_DWARF_MERGER_FRACTION_BY_MASS`), from
the magnetism, the velocity and the rotation of the 204 white dwarfs of the 100 pc sample. Below 0.9 [R5354] reads a second peak of
the mass distribution at 0.80 solar masses, a quarter of the stars of 0.7 to 1 solar mass, and the generator takes that quarter as the
share down to `tuning.WHITE_DWARF_MERGER_FLOOR_SOL`.

The generator's white dwarfs already have the mass the initial-final relation gives (`remnant_mass`), a tail of heavy ones among them, so
a merger product is a white dwarf of that mass that the share flags, not one with a mass of its own: a second draw of the mass would
count the heavy ones twice (a white dwarf of 0.7 solar masses or more is 23 percent of the old ones the relation gives, 11 percent of the
100 pc sample of [R5354]). Over an old population the share comes to 8 to 12 percent of the white dwarfs, inside [R5354]'s 14 percent
with its factor of 2.
"""

from planetgen import tuning


def merger_fraction(mass_sol):
    """The share of the white dwarfs of `mass_sol` solar masses that are double white dwarf merger products, 0 to 1."""
    if mass_sol < tuning.WHITE_DWARF_MERGER_FLOOR_SOL:
        return 0.0
    for low, high, share in tuning.WHITE_DWARF_MERGER_FRACTION_BY_MASS:
        if low <= mass_sol < high:
            return share
    if mass_sol < tuning.WHITE_DWARF_MERGER_FRACTION_BY_MASS[0][0]:
        return tuning.WHITE_DWARF_MERGER_BELOW_SHARE
    return tuning.WHITE_DWARF_MERGER_FRACTION_BY_MASS[-1][2]
