# planetgen/db/compact_pairs.py

"""
The compact binaries of a stored member tree (UX.186, GEN.272), read from what the tree stores: the pair kind from a
node's `interaction` and gas-source bits, the microquasar flag from the accretor's member bits, and the spider kind from
the pulsar's `ablating` bit and its companion's mass. The X-ray luminosity of an accreting pair is worked out on read from
the stored masses, radii (`scatter_members.radius_au`, GEN.952) and orbit by `physics.accretion.x_ray_luminosity_w`, the same call
`generation.multistar.MemberTree.x_ray_luminosities` makes; a row made before the radius was stored has none.
"""

from planetgen.generation import multistar
from planetgen.physics import accretion

KINDS = {
    "high-mass-xray-binary": "High-mass X-ray binary",
    "be-xray-binary": "Be/X-ray binary",
    "low-mass-xray-binary": "Low-mass X-ray binary",
    "black-widow": "Black widow pulsar",
    "redback": "Redback pulsar",
}
"""dict: Each pair kind and its label."""

SPIDER_MASS_LIMIT_SOLAR = 0.1
"""float: A spider's companion under this mass is a black widow, a heavier one a redback (`MemberTree.spider_nodes`)."""


def compact_pairs(members, nodes):
    """
    The compact binaries of a tree, one per node whose two children are single members: an `accreting` pair of a neutron star
    or a black hole and a companion (fed by a wind, a Be star's disc or a donor that fills its lobe) or an ablating pulsar and
    its companion.

    Args:
        members (list[dict]): The tree's members, each with `member_no`, `kind`, `mass_solar`, `parent_node_no`, `flags` and,
            where stored, `radius_au`.
        nodes (list[dict]): The tree's nodes, each with `node_no`, `interaction`, `flags`, `separation_au` and `eccentricity`.

    Returns:
        list[dict]: `node_no`, `kind` (a key of `KINDS`), `compact` and `companion` (the two `member_no`s) and
            `microquasar` (whether the compact object launches jets) and `x_ray_luminosity_w` (watts: an accreting pair's, `None` for a
            spider pair and for a pair whose members have no stored radius), in node order.
    """
    found = []
    for node in sorted(nodes, key=lambda n: n["node_no"]):
        pair = [member for member in members if member["parent_node_no"] == node["node_no"]]
        if len(pair) != 2:
            continue
        spider = next((m for m in pair if (m["flags"] or 0) & multistar.MEMBER_ABLATING), None)
        if spider is not None:
            companion = next(m for m in pair if m is not spider)
            kind = "black-widow" if companion["mass_solar"] < SPIDER_MASS_LIMIT_SOLAR else "redback"
            found.append({"node_no": node["node_no"], "kind": kind, "compact": spider["member_no"],
                          "companion": companion["member_no"], "microquasar": False, "x_ray_luminosity_w": None})
            continue
        compact = next((m for m in pair if m["kind"] in multistar.COMPACT_KINDS), None)
        if node.get("interaction") != "accreting" or compact is None:
            continue
        companion = next(m for m in pair if m is not compact)
        bits = node.get("flags") or 0
        kind = ("be-xray-binary" if bits & multistar.NODE_DISC_FED
                else "high-mass-xray-binary" if bits & multistar.NODE_WIND_FED else "low-mass-xray-binary")
        found.append({"node_no": node["node_no"], "kind": kind, "compact": compact["member_no"],
                      "companion": companion["member_no"],
                      "microquasar": bool((compact["flags"] or 0) & multistar.MEMBER_MICROQUASAR),
                      "x_ray_luminosity_w": _x_ray_luminosity_w(node, compact, companion)})
    return found


def _x_ray_luminosity_w(node, accretor, donor):
    """The pair's X-ray luminosity in watts, or None when either member has no stored radius or the node has no orbit."""
    if accretor.get("radius_au") is None or donor.get("radius_au") is None or not node.get("separation_au"):
        return None
    bits = node.get("flags") or 0
    mode = "wind" if bits & multistar.NODE_WIND_FED else "disc" if bits & multistar.NODE_DISC_FED else "overflow"
    return accretion.x_ray_luminosity_w(
        mode, accretor["kind"], accretor["mass_solar"], accretor["radius_au"], donor["kind"], donor["mass_solar"],
        donor["radius_au"], node["separation_au"], node.get("eccentricity") or 0.0)
