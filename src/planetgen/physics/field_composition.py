# planetgen/physics/field_composition.py

"""
What an interstellar asteroid field made by a collision of two rogue planets is made of (GEN.1102; docs/design/systems/interstellar-asteroid-fields.md sections 3.2 and 7,
docs/design/research/interstellar-field-sources.md point 5)
==========================================================================================================================================================================

Pure functions of two differentiated parents and a collision, no draw. The old `small_body_composition.draw_field("collisional")` picks one taxon at random; this follows the parents.

* **Layers.** A rogue is an ice shell (water-rich ones only), a silicate mantle and an iron-nickel core, from the surface in. `layer_masses` splits the mass; the core fraction is a third
  unless the caller says (IF8).
* **What leaves.** The pieces that fly free are the outermost ones: each parent gives up its own fraction of mass, the unbound share of `rogue_collision.unbound_fraction`, from the
  surface in, ice first, then mantle, then core (`escaping_layers`). The cores are two thirds of the way down, so at the 11 to 14 percent that leaves, a field is ice and silicate
  and holds no iron (`CORE_LEAK_SHARE` 0, row IF23). The cores go into the largest remnant or into the bound debris round it (`remnant_core_fraction`, the average of the two idealised models
  of [R2866] through [R4226], IF8).
* **A collision.** `collision_field` runs the pair through the law of `rogue_collision` (speed to the energy ratio, to the largest remnant, to the unbound share) and gives the field, the
  remnant planet that stays as its own rogue (Boss, 4:54 PM EDT 2026-10-30) and the bound debris. Nothing is a field when the largest remnant keeps over half of the mass: a merger, a hit-and-run or an erosion leaves a remnant and
  what falls back, and `collision_field` gives `None` for the field.
* **Makeup.** `field_makeup` writes the shards as rows of the shared materials table (forsterite and enstatite for the mantle, as the ring does, `kamacite` for iron, `water_ice`) and builds
  a `SmallBodyMakeup` of the kind `asteroid-field` at a porosity the caller draws.

Stand-ins are IF23 to IF25 of `docs/design/research/gap-filling-plan.md` section 6.107. Nothing reads the module yet: the field generation that replaces `draw_field("collisional")`
(GEN.1103, GEN.1104) wires it, and `draw_field("collisional")` stays until then.
"""

import collections
import math

from planetgen.physics import rogue_collision, small_body_composition

DEFAULT_CORE_FRACTION = 1.0 / 3.0
"""float: The core share of the mass of a differentiated rocky rogue (IF8: 0.25 to 0.4, Earth-like 0.33). Tier T1, L1."""
CORE_LEAK_SHARE = 0.0
"""float: The share of the escaping mass that comes from the core while the core is still held under the mantle: none, because the pieces that fly free are the outermost ones. IF23, tier T5,
L0; a hydrocode run of the pair replaces it."""
MANTLE_ROWS = {"forsterite": 0.5, "enstatite": 0.5}
"""dict: The rows of the shared materials table the silicate mantle is written as, mass fractions (peridotite as olivine and pyroxene, the split the ring builder uses). IF24, tier T6, L1."""
CORE_ROW = "kamacite"
"""str: The row the iron-nickel core is written as (4 to 7 percent nickel). IF24."""
ICE_ROW = "water_ice"
"""str: The row the ice shell is written as."""
FIELD_KIND = "asteroid-field"
"""str: The kind of the makeup."""
FIELD_REGIMES = ("catastrophic", "super_catastrophic")
"""tuple: The regimes that leave a field: the largest remnant under a half of the mass, which is the rubble window of the design (22.6 to 30.3 km/s for two Earths and beyond)."""

Parent = collections.namedtuple("Parent", "mass_kg core_fraction ice_fraction")
"""namedtuple: A differentiated rogue: its mass, the share of it in the core and the share in the ice shell (the rest is mantle)."""
Layers = collections.namedtuple("Layers", "ice mantle core")
"""namedtuple: Masses (or shares) of the ice shell, the mantle and the core."""
CollisionOutcome = collections.namedtuple("CollisionOutcome", "regime largest_fraction field_mass_kg field remnant bound_debris_kg")
"""namedtuple: What a collision leaves: the regime, `Mlr / Mtot`, the mass of the field, its `Layers` (masses, `None` when no field), the remnant planet as a `Parent` and the bound debris mass."""


class FieldCompositionError(ValueError):
    """A mass, a fraction or a speed is out of range."""


def _check_parent(parent):
    if not parent.mass_kg > 0.0:
        raise FieldCompositionError(f"the mass {parent.mass_kg!r} is not positive")
    if not (0.0 <= parent.core_fraction < 1.0 and 0.0 <= parent.ice_fraction < 1.0 and parent.core_fraction + parent.ice_fraction < 1.0):
        raise FieldCompositionError("the core and ice fractions must be in 0 to 1 and leave a mantle")


def layer_masses(parent):
    """The `Layers` masses of a parent, kg: its ice shell, its mantle and its core."""
    _check_parent(parent)
    ice = parent.mass_kg * parent.ice_fraction
    core = parent.mass_kg * parent.core_fraction
    return Layers(ice, parent.mass_kg - ice - core, core)


def escaping_layers(parent, loss_fraction, core_leak_share=CORE_LEAK_SHARE):
    """
    The `Layers` masses a parent loses to the field when it gives up `loss_fraction` of its mass from the surface in: the ice shell first, then the mantle, then the core. While the loss is under
    the ice and mantle, `core_leak_share` of it comes from the core instead of the mantle (0 by default). A loss of 1 is the whole body.
    """
    if not 0.0 <= loss_fraction <= 1.0 or not 0.0 <= core_leak_share <= 1.0:
        raise FieldCompositionError("the loss fraction and the leak share are in 0 to 1")
    layers = layer_masses(parent)
    loss = loss_fraction * parent.mass_kg
    core_taken = min(loss * core_leak_share, layers.core)
    loss -= core_taken
    ice_taken = min(loss, layers.ice)
    mantle_taken = min(loss - ice_taken, layers.mantle)
    core_taken += loss - ice_taken - mantle_taken
    return Layers(ice_taken, mantle_taken, min(core_taken, layers.core))


def field_layers(parent1, parent2, unbound_fraction, core_leak_share=CORE_LEAK_SHARE):
    """The `Layers` masses of the field when each parent gives up `unbound_fraction` of its own mass from the surface in: the sum of `escaping_layers` of the two."""
    first = escaping_layers(parent1, unbound_fraction, core_leak_share)
    second = escaping_layers(parent2, unbound_fraction, core_leak_share)
    return Layers(first.ice + second.ice, first.mantle + second.mantle, first.core + second.core)


def remnant_core_fraction(parent1, parent2, largest_fraction):
    """
    The core share of the largest remnant, the average of the two idealised models of IF8: the cores all go into the remnant (as far as it holds them), or the remnant keeps the cores in the
    proportion of the pair. `largest_fraction` is `Mlr / Mtot`; the result is in 0 to 1.
    """
    if not 0.0 < largest_fraction <= 1.0:
        raise FieldCompositionError(f"the largest remnant fraction {largest_fraction!r} is not in (0, 1]")
    total = parent1.mass_kg + parent2.mass_kg
    cores = layer_masses(parent1).core + layer_masses(parent2).core
    remnant = largest_fraction * total
    all_cores = min(cores, remnant) / remnant
    proportional = cores / total
    return 0.5 * (all_cores + proportional)


def _mass_weighted(parent1, parent2, attribute):
    return (getattr(parent1, attribute) * parent1.mass_kg + getattr(parent2, attribute) * parent2.mass_kg) / (parent1.mass_kg + parent2.mass_kg)


def collision_field(parent1, parent2, speed_ms, impact_parameter=0.0, core_leak_share=CORE_LEAK_SHARE):
    """
    The `CollisionOutcome` of a collision of two parents at `speed_ms` and the impact parameter `impact_parameter` (0 head on, 1 grazing, `rogue_collision.critical_impact_parameter` between):
    the energy ratio of `rogue_collision.oblique_disruption_energy_j_kg`, the largest remnant, the regime and, in the catastrophic and super-catastrophic regimes (`FIELD_REGIMES`), the unbound share that makes the field.
    The remnant planet carries `remnant_core_fraction` and the ice of the pair that is left; the debris that is not a field is the mass that neither the remnant nor the field holds.
    """
    if speed_ms <= 0.0 or not 0.0 <= impact_parameter <= 1.0:
        raise FieldCompositionError("the speed must be positive and the impact parameter in 0 to 1")
    _check_parent(parent1)
    _check_parent(parent2)
    target, projectile = (parent1, parent2) if parent1.mass_kg >= parent2.mass_kg else (parent2, parent1)
    total = target.mass_kg + projectile.mass_kg
    radius_t, radius_p = (rogue_collision.combined_radius_m(m.mass_kg) for m in (target, projectile))
    q_r = rogue_collision.specific_impact_energy_j_kg(target.mass_kg, projectile.mass_kg, speed_ms)
    q_star = rogue_collision.oblique_disruption_energy_j_kg(target.mass_kg, projectile.mass_kg, radius_t, radius_p, impact_parameter)
    ratio = q_r / q_star
    largest = rogue_collision.largest_remnant_fraction(ratio)
    grazing = impact_parameter >= rogue_collision.critical_impact_parameter(radius_t, radius_p)
    name = rogue_collision.regime(ratio, grazing=grazing)
    if name not in FIELD_REGIMES:
        remnant = Parent(largest * total, remnant_core_fraction(target, projectile, largest), _mass_weighted(target, projectile, "ice_fraction"))
        return CollisionOutcome(name, largest, 0.0, None, remnant, total - remnant.mass_kg)
    unbound = rogue_collision.unbound_fraction(largest)
    layers = field_layers(target, projectile, unbound, core_leak_share)
    field_mass = layers.ice + layers.mantle + layers.core
    remnant = Parent(largest * total, remnant_core_fraction(target, projectile, largest), _mass_weighted(target, projectile, "ice_fraction"))
    return CollisionOutcome(name, largest, field_mass, layers, remnant, total - remnant.mass_kg - field_mass)


def field_components(layers):
    """The rows of the shared materials table that a field of the masses `layers` is made of, as mass fractions that add up to 1: the ice, the mantle rows and the core row (absent when zero)."""
    total = layers.ice + layers.mantle + layers.core
    if not total > 0.0:
        raise FieldCompositionError("a field with no mass has no composition")
    components = {}
    if layers.ice > 0.0:
        components[ICE_ROW] = layers.ice / total
    for row, share in MANTLE_ROWS.items():
        if layers.mantle > 0.0:
            components[row] = layers.mantle * share / total
    if layers.core > 0.0:
        components[CORE_ROW] = layers.core / total
    return components


def field_makeup(layers, porosity):
    """The `SmallBodyMakeup` of a field of the masses `layers` at the porosity `porosity` (drawn by the caller, `small_body_composition.POROSITY_RANGES['rocky']` for rubble)."""
    if not 0.0 <= porosity < 1.0:
        raise FieldCompositionError(f"the porosity {porosity!r} is not in [0, 1)")
    components = field_components(layers)
    total = layers.ice + layers.mantle + layers.core
    detail = {"source": "collision", "ice_share": layers.ice / total, "mantle_share": layers.mantle / total, "core_share": layers.core / total}
    return small_body_composition.makeup_from_components(FIELD_KIND, components, porosity, detail=detail)
