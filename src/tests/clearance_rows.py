# tests/clearance_rows.py

"""
The rows `physics/clearance.audit` reads, made from a drawn tree (`multistar.draw_member_tree`) or a built system, so a test can ask the
independent verifier whether the stars of a tree clear their orbits without a database (GEN.1015).
"""

import dataclasses

from planetgen.physics import clearance, constants

AU_KM = constants.AU_TO_KM
SOLAR_KG = constants.SOLAR_MASS_TO_KG
HARD_KINDS = ("overlap", "roche", "stability")
"""tuple: The violations that make a tree of stars wrong rather than advisory."""


def rows_of_tree(tree):
    """The star and node rows of a drawn tree, every member a row."""
    stars = [{"id": m.member_no, "mass_kg": m.mass_solar * SOLAR_KG, "radius_km": m.radius_au * AU_KM,
              "parent_node_no": m.parent_node_no, "slot": m.slot} for m in tree.members]
    return {"stars": stars, "nodes": [dataclasses.asdict(node) for node in tree.nodes]}


def rows_of_system(system):
    """The star, dark member and node rows of a built `StarSystem` (stars only, no planets), as the build stored them."""
    slots = system.star_slots or tuple((system.orbit_nodes[0].node_no, slot) for slot in (0, 1))
    stars = [{"id": number, "mass_kg": star.mass, "radius_km": star.radius, "parent_node_no": parent, "slot": slot}
             for number, (star, (parent, slot)) in enumerate(zip(system.stars, slots))]
    stars += [{"id": -1 - number, "mass_kg": dark.mass_solar * SOLAR_KG, "radius_km": (dark.radius_au or 0.0) * AU_KM,
               "parent_node_no": dark.parent_node_no, "slot": dark.slot, "dark": True, "kind": dark.kind} for number, dark in enumerate(system.dark_members)]
    return {"stars": stars, "nodes": [dataclasses.asdict(node) for node in system.orbit_nodes]}


def hard_violations(rows):
    """What the verifier finds wrong in `rows`: overlap, a Roche lobe filled where it should not be, or an unstable triple."""
    return [item for item in clearance.audit(rows) if item.hard and item.kind in HARD_KINDS]
