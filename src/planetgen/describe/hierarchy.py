# planetgen/describe/hierarchy.py

"""
A multiple system's make-up in words, read off its tree of orbits, one short sentence for each pair: "Triple star system. Abralinar A (K9V)
and Abralinar B (K9V) are a close pair. Abralinar C (M1V) orbits Abralinar A and Abralinar B in a wide orbit." The pair's own `Type`
("Binary (K9V/K9V)") says what the pair is; this says what the system is.
"""

COUNT_NAMES = {2: "Binary", 3: "Triple", 4: "Quadruple", 5: "Quintuple", 6: "Sextuple", 7: "Septuple"}
"""dict: How many stars -> what the system is called."""


def _first_mention(star, seen):
    """A star's name, with its spectral type beside it the first time the text names it."""
    if star["name"] in seen:
        return star["name"]
    seen.add(star["name"])
    spectral = (star.get("type") or "").split(" ")[0]
    return f"{star['name']} ({spectral})" if spectral else star["name"]


def _names(items):
    """The names of stars as a list in words: "A", "A and B", "A, B and C"."""
    if len(items) == 1:
        return items[0]
    return ", ".join(items[:-1]) + " and " + items[-1]


def _members(node, stars, nodes):
    """Every star under a node, side by side in the order of the tree (the first side, then the second)."""
    sides = [(star["slot"], [star]) for star in stars if star["node_no"] == node["node_no"]]
    sides += [(other["slot"] or 0, _members(other, stars, nodes)) for other in nodes if other["parent_node_no"] == node["node_no"]]
    return [star for _slot, group in sorted(sides, key=lambda side: side[0]) for star in group]


def _node_sentences(node, stars, nodes, seen):
    """The sentences for one node of the tree and everything under it, or `None` for a node that is not a pair (a young group)."""
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
    sentences = []
    for kind, item in (sides[0], sides[1]):
        if kind == "node":
            inner = _node_sentences(item, stars, nodes, seen)
            if inner is None:
                return None
            sentences += inner
    if sides[0][0] == "star" and sides[1][0] == "star":
        return sentences + [f"{_first_mention(sides[0][1], seen)} and {_first_mention(sides[1][1], seen)} are a {node['kind']} pair."]
    # The second side orbits the first.
    orbiter = [sides[1][1]] if sides[1][0] == "star" else _members(sides[1][1], stars, nodes)
    centre = [sides[0][1]] if sides[0][0] == "star" else _members(sides[0][1], stars, nodes)
    subject = _names([_first_mention(star, seen) for star in orbiter])
    verb = "orbits" if len(orbiter) == 1 else "orbit"
    return sentences + [f"{subject} {verb} {_names([star['name'] for star in centre])} in a {node['kind']} orbit."]


def hierarchy_sentences(stars, nodes):
    """
    The hierarchy of a system of two or more stars as a list of short sentences (the first says what the system is), or `None` for a single
    star or a tree that is not a nest of pairs.

    Args:
        stars (list of dict): One per star, `name`, `type`, `node_no` (the node it hangs from) and `slot`.
        nodes (list of dict): One per orbit node, `node_no`, `parent_node_no`, `slot` and `kind`.
    """
    if len(stars) < 2:
        return None
    root = next((node for node in nodes if node["parent_node_no"] is None), None)
    sentences = _node_sentences(root, stars, nodes, set()) if root else None
    if sentences is None:
        return None
    return [f"{COUNT_NAMES.get(len(stars), f'{len(stars)}-star')} star system."] + sentences
