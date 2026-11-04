# generation/field_topup_build.py

"""
The light asteroid fields of a sector's top-up, built (GEN.1120, GEN.1126; docs/design/systems/region-object-store.md section 8.5).

What the Rogues level leaves unmatched becomes thin clouds of small rubble (`physics/topup_field.plan`): each is an ordinary `AsteroidField` (its composition
and class as any field's) given the planned mass, rms radius, place and velocity (its density is the composition's), and an apparent origin that is not stored as a claim: no
collision is recorded, so the row has no parent.

    made = build(topup.fields, system_config, "region-key:5", address, edge_pc)   # [ToppedField(plan, field, landing)]
"""

import collections

from planetgen.generation import rogue_landing, rogue_topup_build
from planetgen.generation.phenomena.asteroid_field import AsteroidField
from planetgen.physics import constants
from planetgen.util import draw

ToppedField = collections.namedtuple("ToppedField", ("plan", "field", "landing"))
"""One top-up field: the plan's `TopUpField`, the `AsteroidField` and where it is (`rogue_landing.Landing`)."""

def build(fields, system_config, stream_key, address, edge_pc):
    """
    Every field of a top-up, built and placed.

    Args:
        fields (iterable): The plan's `TopUpField`s.
        system_config (SystemConfig): The shared settings.
        stream_key (str): A key unique to this sector's top-up, which the draws are seeded from.
        address (tuple): The `(shell, band, wedge)` of the sector.
        edge_pc (float): The grid's sector edge, pc.

    Returns:
        list: A `ToppedField` for each plan field, in order.
    """
    circular_speed = rogue_topup_build.circular_speed_at(address, edge_pc)
    made = []
    for index, planned in enumerate(fields):
        with draw.bound(f"topup-field:{stream_key}:{index}"):
            field = AsteroidField(system_config)
        field.radius_ly = planned.radius_m / constants.LIGHTYEAR_M
        field.mass_kg = planned.mass_kg
        where = rogue_topup_build.place(planned, address, edge_pc, circular_speed)
        rogue_landing.apply_motion(field, where)
        made.append(ToppedField(planned, field, where))
    return made
