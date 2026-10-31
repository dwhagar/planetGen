# planetgen/physics/mass_boundaries.py

"""
The Mass Boundaries (GEN.933)
=============================

The masses at which one kind of object gives way to another, in one table with a source and a tier for every row. The mass function,
the companion code, the rogue bodies, the star classes and the remnants all read their edges here; before this each typed its own
(0.08 solar masses in the mass function and the companion code against 80 Jupiter masses for the rogue brown dwarfs; 13 Jupiter masses
typed in four modules; the 8 solar mass edge of the core collapse in six).

A row is `Boundary(name, value, unit, source, tier, note)`; `sol(name)` and `kg(name)` give it in solar masses or kilograms. The tier
is the stand-in registry's (`docs/design/research/stand-in-values.md`): a verified figure, or a recalled one until a lane reads the
source. `note` records where the table's rows disagree with one another and which item settles it.

Not here: the display ladder's 0.075 solar mass rung (`util.ladders.MASS`), where masses are shown in solar masses and not in Jupiter
masses, which is Boss's decision (UX.22, 2026-10-10) and stays as he made it. The kinds of remnant a progenitor leaves and their masses
(the neutron star and black hole mass draws, the white dwarf initial-final relation) are GEN.930's, and the hydrogen-burning limit as a
function of [Fe/H] is GEN.936's.
"""

from dataclasses import dataclass

from planetgen.physics import constants, unit_registry as units

_kg_to_mjup = units.converter("kg", "Mjup")
_kg_to_msun = units.converter("kg", "Msun")
_mjup_to_kg = units.converter("Mjup", "kg")
_msun_to_kg = units.converter("Msun", "kg")


@dataclass(frozen=True)
class Boundary:
    """One edge: `value` in `unit` (`Msun` or `Mjup`), where it comes from and how sure the source is."""
    name: str
    value: float
    unit: str
    source: str
    tier: str
    note: str = ""


STAR_BROWN_DWARF_SPLIT_SOL = 0.08
"""float: Where the stellar mass function ends and the brown dwarf one begins (Kroupa 2001's floor of the stars). Both read this one number,
so no mass between them is left out and none is counted twice (GEN.936)."""

BOUNDARIES = (
    Boundary("hydrogen_burning_limit", STAR_BROWN_DWARF_SPLIT_SOL, "Msun",
             "the least mass of a star the mass function draws, and the mass the companion code splits its stars and brown dwarfs at: "
             "the draw split, not the class of the object it draws",
             "Confident Unverified (recalled: Chabrier and Baraffe 2000 give 0.072 solar masses at solar metallicity; "
             "the 0.08 is the round figure the mass function uses)",
             "whether a drawn object burns hydrogen is decided at its place, by its birth [Fe/H] "
             "(`stellar_abundances.hydrogen_burning_limit_sol`: 0.074 at solar Z, 0.092 at Z = 0)."),
    Boundary("deuterium_burning_limit", 13.0, "Mjup",
             "the edge between a planet and a brown dwarf: the working definition of the IAU Working Group on Extrasolar Planets",
             "Confident Unverified (recalled: Spiegel, Burrows and Milsom 2011 give 11 to 16 Jupiter masses with composition)"),
    Boundary("brown_dwarf_upper_limit", STAR_BROWN_DWARF_SPLIT_SOL, "Msun",
             "the heaviest brown dwarf the generator draws: the floor of the stars' mass function, so the two functions meet",
             "Confident Unverified (the draw split above; it was a round 80 Jupiter masses, 0.0764 solar masses, and left the stretch to the "
             "stars' 0.08 empty)",
             "the objects of this stretch burn hydrogen or not by their place's [Fe/H]: the brown dwarf tracks give a stable "
             "late-M or L dwarf to one that does."),
    Boundary("core_collapse_limit", 8.0, "Msun",
             "the initial mass above which a star ends its life in a core collapse (a neutron star or a black hole) and below which it "
             "ends as a white dwarf; also where the giant branch gives way to the supergiants",
             "Confident Unverified (recalled: the textbook 8 solar masses; the real edge depends on metallicity and rotation, 7 to 10)"),
    Boundary("white_dwarf_progenitor_floor", 0.8, "Msun",
             "the lightest star that has become a white dwarf within the age of the Galaxy's disk",
             "Confident Unverified (recalled: a 0.8 solar mass star leaves the main sequence at about 14 Gy)"),
    Boundary("chandrasekhar_limit", constants.CHANDRASEKHAR_LIMIT_SOL, "Msun",
             "the heaviest white dwarf (Chandrasekhar 1931; 1.44 solar masses for a carbon-oxygen white dwarf)",
             "Verified"),
)

SPECTRAL_CLASS_LOWER_MASS_SOL = (("O", 16.0), ("B", 2.1), ("A", 1.4), ("F", 1.04), ("G", 0.8), ("K", 0.45), ("M", 0.0))
"""tuple: The initial mass at which each spectral class begins, heaviest first (the main-sequence masses of the classes' temperatures;
the `star-count` table's rows are by these classes)."""

B_SUBCLASS_LOWER_MASS_SOL = (("B-early", 9.0), ("B-mid", 5.0))
"""tuple: The initial mass at which the B-star sub-classes begin, heaviest first (Table 13's bins); a B star below the last is late B
and keeps the letter's row. `hosts.STAR_COUNT_PROBABILITIES` for mid and early B stars are rederived from these (GEN.449)."""

_BY_NAME = {boundary.name: boundary for boundary in BOUNDARIES}


def boundary(name):
    """The `Boundary` called `name`."""
    return _BY_NAME[name]


def sol(name):
    """The boundary's mass in solar masses."""
    found = _BY_NAME[name]
    if found.unit == "Msun":
        return found.value
    return _kg_to_msun(_mjup_to_kg(found.value))


def jupiter(name):
    """The boundary's mass in Jupiter masses."""
    found = _BY_NAME[name]
    if found.unit == "Mjup":
        return found.value
    return _kg_to_mjup(_msun_to_kg(found.value))


def kg(name):
    """The boundary's mass in kilograms."""
    found = _BY_NAME[name]
    return found.value * (constants.SOLAR_MASS_TO_KG if found.unit == "Msun" else constants.JUPITER_MASS_TO_KG)
