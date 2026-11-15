# tests/test_oxygen_edge_clear_of_cap.py

"""
GEN.1582: the upper oxygen edge of the classes that need a living world's air (M, O, H, P and R) stays clear of the biotic oxygen cap, so a change to either number cannot
drop those classes silently. The cap is `mixture_oxygen.BIOTIC_O2_PRESSURE_MAX_PA`; the edges are the `o2_partial_pa` upper bounds in `physics/data/planet_class_fit.json`.
"""

from planetgen.physics import class_fit, mixture_oxygen

CLASSES_UNDER_THE_CAP = ("M", "O", "H", "P", "R")
EDGE_MARGIN = 0.001
"""float: The least share of the cap by which the edge sits above it (0.1 percent, about 35 Pa: the edge is the cap rounded up to the next 100 Pa)."""


def _oxygen_edges(letter):
    return [c["lt"] for cls in class_fit.CLASSES if cls.letter == letter and cls.kind == "t" for c in cls.when
            if c["quantity"] == "o2_partial_pa" and "lt" in c]


def test_each_class_that_lives_under_the_cap_has_an_oxygen_edge_above_it():
    for letter in CLASSES_UNDER_THE_CAP:
        edges = _oxygen_edges(letter)
        assert edges, f"class {letter} has no upper oxygen edge to hold against the biotic cap"
        for edge in edges:
            assert edge >= mixture_oxygen.BIOTIC_O2_PRESSURE_MAX_PA * (1 + EDGE_MARGIN), (
                f"class {letter}: the upper oxygen edge {edge} Pa is not above the biotic cap {mixture_oxygen.BIOTIC_O2_PRESSURE_MAX_PA} Pa by {EDGE_MARGIN:.1%}")


def test_the_edge_is_not_far_above_the_cap_either():
    # an edge far above the cap would make the class table accept oxygen the biosphere cannot make, and hide a stale number
    for letter in CLASSES_UNDER_THE_CAP:
        for edge in _oxygen_edges(letter):
            assert edge <= mixture_oxygen.BIOTIC_O2_PRESSURE_MAX_PA * 1.01
