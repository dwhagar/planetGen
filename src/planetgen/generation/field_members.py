# planetgen/generation/field_members.py

"""
The stored members of an interstellar asteroid field (GEN.1105, docs/design/systems/interstellar-asteroid-fields.md section 6)
==============================================================================================================================

A field that has a ledger (`AsteroidField.ledger`, the `belt_ledger.Ledger` its mass was solved with: a field the supernova host rule or a collision made) keeps the bodies the ledger stores as
individual objects, heaviest first, by the belt cut rule (every body above `belt_ledger.L1_MASS_KG`, and more of the next heaviest until 90 percent of the field's mass is held, at most 300). A body's rank
is its number, so `<field ID>.<rank hex>` (`object_uid.field_designation`) always names the same body. A field has no star, so a body has no orbit of its own: where it is at any age follows from the
field's seed and `field_dispersal.member_position_m`.

A body takes the field's makeup (the shares of the field's `makeup`: rock, metal, water, ice) and its bulk density (the makeup's, `makeup.bulk_density_kg_m3`, or its grains' times the filling the porosity
leaves), the radius that mass has at that density, and the status word of its roundness (`dwarf planet` where it is round, else `minor body`: the field's own mass clears nothing). A field with no ledger
(a hand-made one, or a Rogues top-up field, whose largest body is under the cut) has no members.

    rows = member_rows(field)       # [{"rank": 1, "mass_kg": ..., "radius_km": ..., "density_g_cm3": ..., "rock_fraction": ..., ..., "status": "minor body"}, ...]
"""

import math

from planetgen.physics import body_status, small_body_composition as sbc


def bulk_density_kg_m3(makeup):
    """The bulk density, kg/m^3, a body of the field has: the makeup's, or its grains' with the void share of the porosity taken out. [Stand-in FD1, tier T9, tweak level L1.]"""
    if makeup.bulk_density_kg_m3 is not None:
        return makeup.bulk_density_kg_m3
    return makeup.grain_density_kg_m3 * (1.0 - makeup.porosity)


def member_rows(field):
    """The rows of the stored members of `field`, heaviest first (none for a field with no ledger or no makeup): `rank`, `mass_kg`, `radius_km`, `density_g_cm3`, the four shares and `status`."""
    ledger = getattr(field, "ledger", None)
    makeup = getattr(field, "makeup", None)
    if ledger is None or makeup is None:
        return []
    density = bulk_density_kg_m3(makeup)
    if not density > 0.0:
        return []
    ice = makeup.water + makeup.ice
    rows = []
    for rank, mass in enumerate(ledger.masses_kg, start=1):
        radius_km = (3.0 * mass / (4.0 * math.pi * density)) ** (1.0 / 3.0) / 1.0e3
        round_body = body_status.roundness(mass, ice) >= 1.0
        rows.append({"rank": rank, "mass_kg": mass, "radius_km": radius_km, "density_g_cm3": density / 1000.0, "rock_fraction": makeup.rock, "metal_fraction": makeup.metal,
                     "water_fraction": makeup.water, "ice_fraction": makeup.ice, "status": body_status.DWARF_PLANET if round_body else body_status.MINOR_BODY})
    return rows
