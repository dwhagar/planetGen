"""A multiple system's make-up in words (describe/hierarchy.py), one short sentence for each pair (UX.253)."""

from planetgen.describe.hierarchy import hierarchy_sentences

STARS = [{"name": "Abralinar A", "type": "K9V Orange Main Sequence Star", "node_no": 0, "slot": 0},
         {"name": "Abralinar B", "type": "K9V Orange Main Sequence Star", "node_no": 0, "slot": 1},
         {"name": "Abralinar C", "type": "M1V Red Main Sequence Star", "node_no": 1, "slot": 1}]
NODES = [{"node_no": 1, "parent_node_no": None, "slot": None, "kind": "wide"},
         {"node_no": 0, "parent_node_no": 1, "slot": 0, "kind": "close"}]


def test_a_triple_says_which_pair_is_close_and_which_star_is_far():
    assert hierarchy_sentences(STARS, NODES) == [
        "Triple star system.",
        "Abralinar A (K9V) and Abralinar B (K9V) are a close pair.",
        "Abralinar C (M1V) orbits Abralinar A and Abralinar B in a wide orbit."]


def test_a_binary_and_a_single_star():
    pair = [dict(NODES[1], parent_node_no=None)]
    assert hierarchy_sentences(STARS[:2], pair)[:2] == ["Binary star system.", "Abralinar A (K9V) and Abralinar B (K9V) are a close pair."]
    assert hierarchy_sentences(STARS[:1], []) is None


def test_two_pairs_read_in_short_sentences():
    stars = [{"name": "Abralinar A", "type": "K9V Orange Main Sequence Star", "node_no": 0, "slot": 0},
             {"name": "Abralinar B", "type": "K9V Orange Main Sequence Star", "node_no": 0, "slot": 1},
             {"name": "Abralinar C", "type": "G2V Yellow Main Sequence Star", "node_no": 2, "slot": 0},
             {"name": "Abralinar D", "type": "G2V Yellow Main Sequence Star", "node_no": 2, "slot": 1}]
    nodes = [{"node_no": 3, "parent_node_no": None, "slot": None, "kind": "wide"},
             {"node_no": 0, "parent_node_no": 3, "slot": 0, "kind": "close"},
             {"node_no": 2, "parent_node_no": 3, "slot": 1, "kind": "close"}]
    assert hierarchy_sentences(stars, nodes) == [
        "Quadruple star system.",
        "Abralinar A (K9V) and Abralinar B (K9V) are a close pair.",
        "Abralinar C (G2V) and Abralinar D (G2V) are a close pair.",
        "Abralinar C and Abralinar D orbit Abralinar A and Abralinar B in a wide orbit."]


def test_a_pair_orbited_by_a_star_that_is_orbited_in_turn_reads_in_order():
    stars = [{"name": "A", "type": "K9V", "node_no": 0, "slot": 0}, {"name": "B", "type": "K9V", "node_no": 0, "slot": 1},
             {"name": "C", "type": "M1V", "node_no": 1, "slot": 1}, {"name": "D", "type": "M5V", "node_no": 2, "slot": 1}]
    nodes = [{"node_no": 2, "parent_node_no": None, "slot": None, "kind": "wide"},
             {"node_no": 1, "parent_node_no": 2, "slot": 0, "kind": "medium"},
             {"node_no": 0, "parent_node_no": 1, "slot": 0, "kind": "close"}]
    assert hierarchy_sentences(stars, nodes) == [
        "Quadruple star system.", "A (K9V) and B (K9V) are a close pair.", "C (M1V) orbits A and B in a medium orbit.",
        "D (M5V) orbits A, B and C in a wide orbit."]
