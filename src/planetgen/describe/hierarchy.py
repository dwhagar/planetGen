# planetgen/describe/hierarchy.py

"""
A multiple system's make-up in words, read off its tree of orbits: "Triple star system: a close pair of Abralinar A (K9V)
and Abralinar B (K9V), orbited by Abralinar C (M1V) in a wide orbit." The pair's own `Type` ("Binary (K9V/K9V)") says what
the pair is; this says what the system is.
"""

COUNT_NAMES = {2: "Binary", 3: "Triple", 4: "Quadruple", 5: "Quintuple", 6: "Sextuple", 7: "Septuple"}
"""dict: How many stars -> what the system is called."""


def _star_text(star):
    spectral = (star.get("type") or "").split(" ")[0]
    return f"{star['name']} ({spectral})" if spectral else star["name"]


def _node_text(node, stars, nodes):
    """One node of the tree in words, or `None` for a node that is not a pair (a young group)."""
    if node["kind"] == "group":
        return None
    sides = {}
    for star in stars:
        if star["node_no"] == node["node_no"]:
            sides[star["slot"]] = ("star", star)
    for other in nodes:
        if other["parent_node_no"] == node["node_no"]:
            sides[other["slot"] or 0] = ("node", other)
    if set(sides) != {0, 1}:
        return None
    parts = [_star_text(item) if kind == "star" else _node_text(item, stars, nodes) for kind, item in (sides[0], sides[1])]
    if None in parts:
        return None
    if sides[0][0] == "star" and sides[1][0] == "star":
        return f"a {node['kind']} pair of {parts[0]} and {parts[1]}"
    if sides[0][0] == "node" and sides[1][0] == "star":
        return f"{parts[0]}, orbited by {parts[1]} in a {node['kind']} orbit"
    if sides[0][0] == "star":
        return f"{parts[0]}, orbited by {parts[1]} in a {node['kind']} orbit"
    return f"{parts[0]}, orbited by {parts[1]}"


def hierarchy_text(stars, nodes):
    """
    The hierarchy of a system of two or more stars as a sentence without its full stop, or `None` for a single star or a
    tree that is not a nest of pairs.

    Args:
        stars (list of dict): One per star, `name`, `type`, `node_no` (the node it hangs from) and `slot`.
        nodes (list of dict): One per orbit node, `node_no`, `parent_node_no`, `slot` and `kind`.
    """
    if len(stars) < 2:
        return None
    root = next((node for node in nodes if node["parent_node_no"] is None), None)
    text = _node_text(root, stars, nodes) if root else None
    if text is None:
        return None
    return f"{COUNT_NAMES.get(len(stars), f'{len(stars)}-star')} star system: {text}"
