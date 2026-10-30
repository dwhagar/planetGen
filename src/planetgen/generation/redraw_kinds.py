# planetgen/generation/redraw_kinds.py

"""
Redraw kinds
============

Every place the generator throws a drawn item away, or moves or cuts it, to make it fit has a named kind (TEST.476, GEN.1006). The
place calls `draw.redrawn(kind)` (which reads no stream, so what is generated is the same counted or not) and the census of
`tests/test_redraw_rates.py` reports each kind's rate per generated system, by number of stars and placed-star class, against the
ceilings in `tests/data/redraw_rate_expectations.json`.

A kind has a `scope` (what is redrawn), an `overlap` flag (the kinds the overlap-free design of GEN.1006 closes: their count is held at
exactly zero once GEN.1015 deletes the GEN.1004 quick fix, so the test holds it as an invariant and not under a ceiling) and a one-line
`summary`. Not counted: the rejection samplers of a sector's star population (`star_population`, `bright_stars`), which accept or reject
a draw by its weight and are not a fit to a neighbour, and the bisections and fixed-point loops that draw nothing.
"""

from collections import namedtuple

Kind = namedtuple("Kind", "scope overlap summary")

KINDS = {
    # --- the pair placement of a star tree (multistar._Builder.place): one count per try that is thrown away ---------------------------
    "pair.axis_out_of_range": Kind("star-tree", False, "an orbit axis outside the stable window of its parent or the place's widest axis"),
    "pair.no_allowed_period": Kind("star-tree", False, "no period the law allows for this pair between the Roche floor and the ceiling (the legal reason)"),
    "pair.overlap": Kind("star-tree", True, "a pair whose closest approach does not clear what it holds (the Roche limit or the inner apocentres)"),
    "pair.interaction_gate": Kind("star-tree", False, "an interaction the pair's members or the age do not allow (a compact accretor, a contact binary of the wrong age)"),
    "pair.unstable_triple": Kind("star-tree", False, "a triple that fails the hierarchical stability screen"),
    "pair.kozai": Kind("star-tree", False, "a triple the Kozai-Lidov screen calls unsafe"),
    "pair.gravitational_waves": Kind("star-tree", False, "a compact or double white dwarf pair that would merge within the age by gravitational waves"),
    "pair.brown_dwarf_desert": Kind("star-tree", False, "a brown dwarf closer than the desert's edge, thinned to the observed share"),
    # --- the whole tree (multistar._grow_tree and _contact_pair_alone) ------------------------------------------------------------------
    "tree.restart_pair_failed": Kind("star-tree", False, "a tree drawn again because one of its pairs found no orbit in the retries"),
    "tree.restart_overlap": Kind("star-tree", True, "a tree drawn again because the final extents check found a pair that does not clear what it holds"),
    "tree.member_dropped": Kind("star-tree", False, "the lightest loose member dropped after every restart failed"),
    "tree.drawn_at_earlier_age": Kind("star-tree", False, "a tree drawn again at the end of the placed star's lifetime because the star had collapsed by its age"),
    "tree.second_star_forced": Kind("star-tree", False, "a tree drawn again because a companion was wanted and the first draw had none"),
    "star.mass_lightened": Kind("star", False, "a member lightened by a fifth because it had already collapsed at the system's age"),
    "accreting.donor": Kind("star-tree", False, "a donor of an X-ray, cataclysmic or AM CVn binary drawn again because it cannot fill its lobe in the period window"),
    # --- the companion set of a system (system.StarSystem.__init__ and _companions_fit) ------------------------------------------------
    "companions.star_count": Kind("system", False, "the companions drawn again because the star count is not the one the directives force"),
    "companions.pair_kind": Kind("system", False, "the companions drawn again because the pair is not the close or wide kind the directives ask for"),
    "companions.overlap": Kind("system", True, "the companions drawn again because a star as built is wider than its orbit (the quick fix that redraws the whole set)"),
    "companions.circumbinary_floor": Kind("system", False, "the companions drawn again because a close pair leaves no habitable zone outside its stability limit"),
    # --- the planets of a system -------------------------------------------------------------------------------------------------------
    "system.planets_regenerated": Kind("system", False, "every planet drawn again because the forced habitable world, belt or class is not there"),
    "system.planets_regenerated_after_age": Kind("system", False, "every planet drawn again once more under the settled age of the system"),
    "planet.request_missed": Kind("planet", False, "a planet drawn again because it is not the class or habitable world a slot asked for"),
    "planet.barred_world": Kind("planet", False, "a planet drawn again because it came out habitable where habitable worlds are barred"),
    "planet.hill_spacing_push": Kind("planet", False, "a planet moved outward to clear the mutual Hill radius of the one inside it"),
    "planet.belt_spacing_push": Kind("planet", False, "a planet moved outward to clear the belt inside it"),
    "belt.spacing_push": Kind("belt", False, "a belt moved outward to clear the planet or belt inside it"),
    "planet.circumbinary_floor_push": Kind("planet", False, "the innermost body moved out past a close pair's stability limit"),
    "planet.orbit_ceiling_trim": Kind("planet", False, "an outer planet or belt cut because the spacing pass pushed it past the orbit ceiling"),
    "planet.cross_star_trim": Kind("planet", False, "an outer planet cut so the two stars of a wide pair do not encroach on each other's planets"),
    "moon.too_large": Kind("moon", False, "a moon drawn again because it would be larger than its parent could hold"),
    "comet.ceiling": Kind("comet", False, "a comet drawn again because its orbit lies past the star's ceiling"),
    "comet.composition": Kind("comet", False, "a comet's CO and CO2 shares drawn again because one lies outside the observed bounds"),
    "giant.heavy_mass": Kind("giant", False, "a giant's heavy-element mass drawn again because it lies outside the allowed fraction of its mass"),
    "snow_line.scatter": Kind("system", False, "a snow-line scatter drawn again because it lies beyond the sigma cap"),
    "star.model_redraw": Kind("star", False, "a primary star drawn again because it breaks a directive (mass, age, remnant, luminosity) or is already dead"),
}
"""dict: Every redraw kind `draw.redrawn` may be called with, by name."""

OVERLAP_KINDS = tuple(sorted(name for name, kind in KINDS.items() if kind.overlap))
"""tuple: The kinds GEN.1006 closes; the test holds their count at exactly zero once GEN.1015 has deleted the quick fix."""
