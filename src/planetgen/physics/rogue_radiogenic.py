# planetgen/physics/rogue_radiogenic.py

"""
The long-lived heat of a thrown-out planet, from the star it left (GEN.1134)
============================================================================

A rogue rock drew its thorium and uranium from a host-less offset of 0.22 dex (`radiogenic.draw_content(None)`), a potassium and aluminium-26 content of its own, and was called
geologically active if its heat flow reached 0.04 W/m2 (`ROGUE_ACTIVE_HEAT_FLUX_W_M2`), a number between Mars's 0.02 and the Earth's 0.087 (`rogue-planet-review.md` row 12). A planet that was
ejected was made of its star's gas, so it carries the same rock-forming ratios a bound planet of that star does. This module reads them the way a bound planet does, from the ejection record:

* **The offset** of thorium and uranium from the solar scale is `radiogenic.host_offset_dex` for the star's [Fe/H], age and the key of its draws (`stellar_composition.abundance_key`): one number for the star,
  the same for every planet and rogue it made.
* **The r-process enrichment** of the star's birth gas is `rprocess_enrichment.th_u_factor_of_shape` for the star's age in the galaxy it was born in, so an old star's rogues carry less thorium over uranium
  than a young star's, as its bound planets do.
* **The content** (`content`) is `radiogenic.draw_content` with those two, the same draws in the same order, so a rogue and a bound rocky planet of one star, given the same random stream, have
  the same concentrations and the same aluminium-26 factor. The rogue has no star lighting it after the ejection but its rock still decays: its heat at its age is `interior_heat.heat_flow` on that content.
* **Geologically active** (`is_active`) stays a heat-flow cut, now held here with its tier. A physical state was tried and does not discriminate: the one-box mantle of `interior_heat` stays above
  `tidal_heating.SOLIDUS_K` (1,700 K) at every age and mass from the Moon's to 3 Earth masses, 1,950 K and up, so "the mantle still melts" is true of every body. The cut is between Mars (about 0.02 W/m2 today,
  recalled) and the Earth (0.087), at 0.04, a hand-set value (tier T6, level L3) that waits for a model of when a stagnant lid stops erupting.

Pure functions; the random numbers come from the caller. Nothing calls them yet (GEN.1111 wires them). Stand-ins are in `docs/design/research/gap-filling-plan.md` section 6.114.
"""

from planetgen.physics import interior_heat, radiogenic, rprocess_enrichment, stellar_composition
from planetgen.util import draw

ACTIVE_HEAT_FLUX_W_M2 = 0.04
"""float: A rocky rogue whose heat flow reaches this is still geologically active: between Mars (about 0.02 W/m2) and the Earth (0.087). Source tier T6, tweak level L3."""


class RogueRadiogenicError(ValueError):
    """An input the model cannot use."""


def host_key(host_mass_msun=None, host_distance_ly=None, uid=None, group=None):
    """The key of the host's composition draws (`stellar_composition.abundance_key`): the same one its bound planets use."""
    return stellar_composition.abundance_key(uid=uid, mass=host_mass_msun, distance_ly=host_distance_ly, group=group)


def host_offset_dex(host_fe_h, host_age_gy, key):
    """The thorium and uranium offset of the host's rock from the solar scale, dex (`radiogenic.host_offset_dex`)."""
    if not host_age_gy > 0.0:
        raise RogueRadiogenicError("the host's age must be positive")
    return radiogenic.host_offset_dex(float(host_fe_h), float(host_age_gy), key)


def th_u_factor(host_age_gy, shape=None):
    """`F` of the host's birth gas for its age in the galaxy of `shape` (the default galaxy's when none is given), at least 1."""
    from planetgen.galaxy import shape_draw
    if not host_age_gy > 0.0:
        raise RogueRadiogenicError("the host's age must be positive")
    return rprocess_enrichment.th_u_factor_of_shape(shape or shape_draw.default_shape(None), round(float(host_age_gy), 6))


def content(host_fe_h, host_age_gy, key, shape=None, refractory_enrichment=1.0, rng=draw):
    """
    The `radiogenic.Content` of a rogue rock made of the host star's gas, at formation.

    Args:
        host_fe_h, host_age_gy (float): The ejecting star's [Fe/H] and age.
        key (str): The key of its composition draws (`host_key`).
        shape: The `GalaxyShape` it was born in (the default galaxy's when `None`).
        refractory_enrichment (float): How many times the Earth's share the rock's refractory elements are (`radiogenic.draw_content`); 1 when its aluminium is not known.
        rng: The draw source.
    """
    return radiogenic.draw_content(host_offset_dex(host_fe_h, host_age_gy, key), rng=rng, refractory_enrichment=refractory_enrichment,
                                   th_u_factor=th_u_factor(host_age_gy, shape))


def heat(content_at_formation, interior, mass_kg, radius_km, age_gy, rock):
    """The rogue's `interior_heat.Heat` at its age: no star warms it, the rock's own decay and the heat it was born with are all there is."""
    return interior_heat.heat_flow(content_at_formation, interior, mass_kg, radius_km, age_gy, rock)


def is_active(flux_w_m2):
    """Whether a rocky body is still geologically active: its heat flow reaches `ACTIVE_HEAT_FLUX_W_M2`."""
    return flux_w_m2 >= ACTIVE_HEAT_FLUX_W_M2
