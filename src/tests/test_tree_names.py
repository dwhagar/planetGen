# tests/test_tree_names.py

"""
GEN.312: the names of a tree. Path codes (`Aa`, `Ab`, `Ba`) inside a branch, a branch word at each child of a wide split, planet stems that carry the host's code
(`Pikkita Ba I`), every member of a generated tree with a unique, deterministic name, and an old two-star system keeping its names.
"""

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.names import tree_names
from planetgen.util import draw


def _node(node_no, parent, slot, kind):
    return {"node_no": node_no, "parent_node_no": parent, "slot": slot, "kind": kind}


# the design's seven stars: a wide root; branch X = ((x1, x2), x3), branch Y = ((y1, y2), (y3, y4))
SEVEN = ([_node(0, None, 0, "wide"), _node(1, 0, 0, "close"), _node(2, 1, 0, "close"), _node(3, 0, 1, "close"), _node(4, 3, 0, "close"), _node(5, 3, 1, "close")],
         [(2, 0), (2, 1), (1, 1), (4, 0), (4, 1), (5, 0), (5, 1)])


def test_the_seven_star_example_of_the_design():
    nodes, members = SEVEN
    assert tree_names.word_slots(nodes, members) == [(0, 0), (0, 1)]
    names = tree_names.name_tree("Voranthis", nodes, members, ("Kelmoor", "Pikkita"))
    assert names.members == ("Voranthis Kelmoor Aa", "Voranthis Kelmoor Ab", "Voranthis Kelmoor B",
                             "Voranthis Pikkita Aa", "Voranthis Pikkita Ab", "Voranthis Pikkita Ba", "Voranthis Pikkita Bb")
    # the branch that holds the primary uses the first word for its planets, any other branch its own word; a coded star adds its code
    assert names.member_stems == ("Voranthis Aa", "Voranthis Ab", "Voranthis B", "Pikkita Aa", "Pikkita Ab", "Pikkita Ba", "Pikkita Bb")
    # the pair about a node: the branch word at the branch root, else the branch word and the code
    assert names.nodes[1] == "Voranthis Kelmoor" and names.nodes[2] == "Voranthis Kelmoor A" and names.nodes[5] == "Voranthis Pikkita B"
    assert names.node_stems[1] == "Voranthis" and names.node_stems[2] == "Voranthis A" and names.node_stems[3] == "Pikkita" and names.node_stems[5] == "Pikkita B"


def test_a_tree_with_no_wide_split_is_one_word_with_codes():
    nodes = [_node(0, None, 0, "close"), _node(1, 0, 1, "close")]
    names = tree_names.name_tree("Voranthis", nodes, [(0, 0), (1, 0), (1, 1)])
    assert names.members == ("Voranthis A", "Voranthis Ba", "Voranthis Bb")
    assert names.member_stems == names.members and names.node_stems == {0: "Voranthis", 1: "Voranthis B"}
    assert tree_names.word_slots(nodes, [(0, 0), (1, 0), (1, 1)]) == []


def test_the_codes_go_letter_letter_number_then_dotted_numbers():
    nodes = [_node(0, None, 0, "close"), _node(1, 0, 0, "close"), _node(2, 1, 0, "close"), _node(3, 2, 0, "close")]
    members = [(3, 0), (3, 1), (2, 1), (1, 1), (0, 1)]
    names = tree_names.name_tree("Voranthis", nodes, members).members
    assert names == ("Voranthis Aa1.1", "Voranthis Aa1.2", "Voranthis Aa2", "Voranthis Ab", "Voranthis B")


def test_a_group_node_goes_on_past_two_children():
    nodes = [_node(0, None, 0, "group")]
    assert tree_names.name_tree("Voranthis", nodes, [(0, n) for n in range(4)]).members == tuple(f"Voranthis {letter}" for letter in "ABCD")


def test_two_stars_are_named_as_they_always_were():
    close = tree_names.name_tree("Voranthis Prime", [_node(0, None, 0, "close")], [(0, 0), (0, 1)])
    assert close.members == ("Voranthis Prime A", "Voranthis Prime B") and close.node_stems == {0: "Voranthis Prime"}
    wide = tree_names.name_tree("Voranthis Prime", [_node(0, None, 0, "wide")], [(0, 0), (0, 1)], ("Green", "Little"))
    assert wide.members == ("Voranthis Green", "Voranthis Little") and wide.member_stems == ("Voranthis", "Little")


def test_a_lone_star_has_the_system_name():
    assert tree_names.name_tree("Voranthis", [], [(None, 0)]).members == ("Voranthis",)


def test_a_dark_member_has_a_code_among_the_stars():
    nodes = [_node(0, None, 0, "close"), _node(1, 0, 0, "close")]
    assert tree_names.name_tree("Voranthis", nodes, [(1, 0), (1, 1), (0, 1)]).members == ("Voranthis Aa", "Voranthis Ab", "Voranthis B")


def test_the_branch_words_are_needed_and_read_back():
    nodes, members = SEVEN
    with pytest.raises(ValueError):
        tree_names.name_tree("Voranthis", nodes, members, ("Kelmoor",))
    names = tree_names.name_tree("Voranthis Prime", nodes, members, ("Kelmoor", "Pikkita")).members
    assert tree_names.words_from_names("Voranthis Prime", nodes, members, list(names)) == ("Kelmoor", "Pikkita")
    assert tree_names.words_from_names("Voranthis Prime", nodes, members, [*names[:6], "Voranthis Pikkita Ba Bb"]) is None
    assert tree_names.words_from_names("Voranthis Prime", nodes, members, [*names[:2], "Voranthis Other B", *names[3:]]) is None
    renamed = ["Voranthis Kelmoor Aa", "Voranthis Kelmoor Ab", "Voranthis Kelmoor B", "Voranthis Pikkita Aa", "Voranthis Pikkita Ab", "Voranthis Pikkita Ba", "Voranthis Kelmoor Bb"]
    assert tree_names.words_from_names("Voranthis Prime", nodes, members, renamed) is None
    assert tree_names.words_from_names("Voranthis", [_node(0, None, 0, "close")], [(0, 0), (0, 1)], ["Voranthis A", "Voranthis B"]) == ()


def _system(seed, stars, wide=None, planets=False):
    config = SystemConfig()
    config.STARS = stars
    config.PLANETS = planets
    if wide is not None:
        config.WIDE_BINARY = wide
    with draw.bound(seed):
        return StarSystem(config)


@pytest.mark.parametrize("stars", [3, 4, 5])
def test_every_member_of_a_generated_tree_has_a_unique_name_and_two_builds_agree(stars):
    checked = 0
    for seed in range(12):
        system = _system(seed, stars)
        again = _system(seed, stars)
        assert [star.name for star in system.stars] == [star.name for star in again.stars]
        assert system.star_words == again.star_words
        names = [star.name for star in system.stars] + list(system.dark_member_names())
        assert len(names) == len(set(names)), names
        nodes, members = system._name_inputs()
        assert len(system.star_words or ()) == len(tree_names.word_slots(nodes, members))
        assert tree_names.words_from_names(system.name, nodes, members, [star.name for star in system.stars] + [None] * len(system.dark_members)) in (
            system.star_words or (), None) or not system.star_words
        checked += len(system.stars) > 2
    assert checked > 6


def test_the_planets_of_a_generated_triple_carry_the_stem_of_their_host():
    named = 0
    for seed in range(20):
        system = _system(seed, 3, planets=True)
        names = [body.name for body in system.all_planets if body.body_type != "a"]
        assert len(names) == len(set(names))
        nodes, members = system._name_inputs()
        tree = tree_names.name_tree(system.name, nodes, members, system.star_words or ())
        stems = set(tree.member_stems) | set(tree.node_stems.values())
        for name in names:
            assert name.rpartition(" ")[0] in stems, (name, stems)
            named += 1
    assert named > 10
