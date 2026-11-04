# planetgen/physics/mass_loss.py

"""
Orbits around a star that lost mass (GEN.318; docs/design/systems/exotic-environments-planets-and-compact-binaries.md 1.2)
===================================================================================================================

A star of M_initial solar masses that becomes a white dwarf of M_final loses the difference on its giant branches, slowly
and in every direction. Two things follow for a planet that was there:

* **Engulfment.** The giant's envelope took, or drove off, every planet inside `engulfment_radius_au` of the progenitor: a lookup through
  `WD_ENGULFMENT_POINTS_AU`, the Jovian-planet points of [R3225], ApJ 761, 121 (2.7 AU for a
  solar-mass star, 10 AU for five solar masses).
* **Expansion.** A planet outside it keeps its angular momentum while the star loses mass adiabatically, so
  `a (M_star + M_p)` is constant ([R1745, R2149, R4534]; section 1.2 of the
  note): the orbit grows by `(M_initial + M_p) / (M_final + M_p)`, 1.85 for a star of one solar mass that leaves 0.54.

The drawn distances of a white dwarf's planets are the ones the planets formed at; this is what happens to them.
"""

from planetgen import tuning
from planetgen.physics import remnant_mass, unit_registry as units

_kg_to_msun = units.converter("kg", "Msun")
_msun_to_kg = units.converter("Msun", "kg")


def progenitor_mass_sol(remnant_mass_sol, initial_mass_sol=None):
    """The mass, solar masses, of the star a white dwarf of `remnant_mass_sol` was born as: `initial_mass_sol` when
    the model gave one, else the initial-final relation inverted (`remnant_mass.progenitor_mass_sol`), kept inside
    `WD_PROGENITOR_MASS_RANGE_SOL`."""
    if initial_mass_sol is not None:
        return initial_mass_sol
    low, high = tuning.WD_PROGENITOR_MASS_RANGE_SOL
    return min(high, max(low, remnant_mass.progenitor_mass_sol(remnant_mass_sol)))


def engulfment_radius_au(progenitor_sol):
    """The distance, AU, inside which a progenitor of `progenitor_sol` solar masses removed every planet: straight lines
    through `WD_ENGULFMENT_POINTS_AU`, the first value below the first point and the last segment's slope above the last."""
    points = tuning.WD_ENGULFMENT_POINTS_AU
    if progenitor_sol <= points[0][0]:
        return points[0][1]
    for (mass_low, au_low), (mass_high, au_high) in zip(points, points[1:]):
        if progenitor_sol <= mass_high:
            return au_low + (au_high - au_low) * (progenitor_sol - mass_low) / (mass_high - mass_low)
    (mass_low, au_low), (mass_high, au_high) = points[-2:]
    return au_high + (au_high - au_low) / (mass_high - mass_low) * (progenitor_sol - mass_high)


def expansion_factor(initial_mass_kg, final_mass_kg, planet_mass_kg=0.0):
    """How much an orbit grows when the star it circles loses mass slowly and evenly: `(M_initial + M_p) / (M_final +
    M_p)` (1 when the star lost none)."""
    return (initial_mass_kg + planet_mass_kg) / (final_mass_kg + planet_mass_kg)


def member_formation_mass_kg(member):
    """The mass, kg, a star had when its planets formed: a white dwarf's progenitor's, any other star's own."""
    if getattr(member, "yerkes_class", None) in ("VII", "D"):
        progenitor = progenitor_mass_sol(_kg_to_msun(member.mass),
                                         getattr(member, "initial_mass_sol", None))
        return max(_msun_to_kg(progenitor), member.mass)
    return member.mass
