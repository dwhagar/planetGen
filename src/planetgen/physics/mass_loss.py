# planetgen/physics/mass_loss.py

"""
Orbits around a star that lost mass (GEN.318; docs/design/systems/exotic-environments-planets-and-compact-binaries.md 1.2)
===================================================================================================================

A star of M_initial solar masses that becomes a white dwarf of M_final loses the difference on its giant branches, slowly
and in every direction. Two things follow for a planet that was there:

* **Engulfment.** The giant's envelope took, or drove off, every planet inside `engulfment_radius_au` of the progenitor:
  `WD_ENGULFMENT_BASE_AU` for a solar-mass star and `WD_ENGULFMENT_AU_PER_PROGENITOR_MSUN` more per solar mass above
  one.
* **Expansion.** A planet outside it keeps its angular momentum while the star loses mass adiabatically, so
  `a (M_star + M_p)` is constant (Jeans 1924; Hadjidemetriou 1963; Veras et al. 2011, MNRAS 417, 2104; section 1.2 of the
  note): the orbit grows by `(M_initial + M_p) / (M_final + M_p)`, 1.85 for a star of one solar mass that leaves 0.54.

The drawn distances of a white dwarf's planets are the ones the planets formed at; this is what happens to them.
"""

from planetgen import tuning
from planetgen.physics import constants


def progenitor_mass_sol(remnant_mass_sol, initial_mass_sol=None):
    """The mass, solar masses, of the star a white dwarf of `remnant_mass_sol` was born as: `initial_mass_sol` when
    the model gave one, else the remnant relation inverted (`WHITE_DWARF_REMNANT_MASS_SOL`), kept inside
    `WD_PROGENITOR_MASS_RANGE_SOL`."""
    if initial_mass_sol is not None:
        return initial_mass_sol
    relation = tuning.WHITE_DWARF_REMNANT_MASS_SOL
    low, high = tuning.WD_PROGENITOR_MASS_RANGE_SOL
    return min(high, max(low, (remnant_mass_sol - relation["intercept"]) / relation["slope"]))


def engulfment_radius_au(progenitor_sol):
    """The distance, AU, inside which a progenitor of `progenitor_sol` solar masses removed every planet."""
    return tuning.WD_ENGULFMENT_BASE_AU + tuning.WD_ENGULFMENT_AU_PER_PROGENITOR_MSUN * (progenitor_sol - 1.0)


def expansion_factor(initial_mass_kg, final_mass_kg, planet_mass_kg=0.0):
    """How much an orbit grows when the star it circles loses mass slowly and evenly: `(M_initial + M_p) / (M_final +
    M_p)` (1 when the star lost none)."""
    return (initial_mass_kg + planet_mass_kg) / (final_mass_kg + planet_mass_kg)


def member_formation_mass_kg(member):
    """The mass, kg, a star had when its planets formed: a white dwarf's progenitor's, any other star's own."""
    if getattr(member, "yerkes_class", None) in ("VII", "D"):
        progenitor = progenitor_mass_sol(member.mass / constants.SOLAR_MASS_TO_KG,
                                         getattr(member, "initial_mass_sol", None))
        return max(progenitor * constants.SOLAR_MASS_TO_KG, member.mass)
    return member.mass
