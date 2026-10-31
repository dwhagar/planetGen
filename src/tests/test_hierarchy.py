"""A multiple system's make-up in words (describe/hierarchy.py)."""

from planetgen.describe.hierarchy import hierarchy_text

STARS = [{"name": "Abralinar A", "type": "K9V Orange Main Sequence Star", "node_no": 0, "slot": 0},
         {"name": "Abralinar B", "type": "K9V Orange Main Sequence Star", "node_no": 0, "slot": 1},
         {"name": "Abralinar C", "type": "M1V Red Main Sequence Star", "node_no": 1, "slot": 1}]
NODES = [{"node_no": 1, "parent_node_no": None, "slot": None, "kind": "wide"},
         {"node_no": 0, "parent_node_no": 1, "slot": 0, "kind": "close"}]


def test_a_triple_says_which_pair_is_close_and_which_star_is_far():
    assert hierarchy_text(STARS, NODES) == (
        "Triple star system: a close pair of Abralinar A (K9V) and Abralinar B (K9V), orbited by Abralinar C (M1V) in a wide orbit")


def test_a_binary_and_a_single_star():
    pair = [dict(NODES[1], parent_node_no=None)]
    assert hierarchy_text(STARS[:2], pair).startswith("Binary star system: a close pair of")
    assert hierarchy_text(STARS[:1], []) is None
