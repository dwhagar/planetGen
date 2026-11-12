# tests/test_materials_reachability.py

"""
TEST.685: every row of the materials table is drawn by some path, or tagged vocabulary-only.

The table (`physics/materials.py`, `physics/data/materials.json`) is the one list of what a body can contain, and a row nobody reads is data that looks like a feature. A row is
*reached* when a module that uses the table, or `tuning.py` (whose component lists the table resolves), names it by id, name or alias, or when a reached group lists it. A row that
is not reached must be listed in `VOCABULARY_ONLY` (built from `scripts/materials_data.py`) with why; an unlisted orphan fails, and so does a listed row that has since been reached,
and a listed name that is no row (the loader refuses that one). The rows a new item adds carry their tag with the item (GEN.1355 to GEN.1357).

The reach is a literal-string search over the modules, not a call graph: a name built in an f-string is missed, and the row must then be tagged or named. The test also lists the
doubles it finds (the same formula under a grain-size or powder variant of the row), and fails on one that is not acknowledged below.
"""

import ast
import pathlib
import re

from planetgen.physics import materials

SOURCE = pathlib.Path(__file__).resolve().parents[1] / "planetgen"
TABLE_MODULE = pathlib.Path("physics") / "materials.py"
ALWAYS_SCANNED = (pathlib.Path("tuning.py"), pathlib.Path("physics/rock_types.py"), pathlib.Path("physics/evaporite_sequence.py"))
"""tuple: Modules scanned for names though they do not use the table: `tuning.ASTEROID_COMPONENTS` is the vocabulary the table resolves, the CIPW norm (`rock_types.py`) names the minerals it makes (nepheline, leucite, hypersthene, wollastonite, periclase), and the evaporite sequence (`evaporite_sequence.py`, GEN.1366) names the salts, carbonates and perchlorates it dries out of a brine; it takes the molar masses from its own formula table, so it does not import the materials table and the scan would otherwise miss it."""

VARIANT_SUFFIX = re.compile(r"_(powder|fine|lt\d+um)$")
"""re.Pattern: How a row id says it is a form of another material (a grain size, a powder)."""

ACKNOWLEDGED_DOUBLES = {
    frozenset({"hematite", "hematite_fine", "hematite_lt10um"}): "hematite in three grain sizes: each carries its own band albedo from a different source",
    frozenset({"anorthite", "anorthosite_powder"}): "the mineral and the powdered rock sample whose measured albedo the ground recipes use",
}
"""dict: The doubles the test knows of, with why each row stays. A new double fails until it is merged away or named here."""


def _uses_table(tree):
    for node in ast.walk(tree):
        if isinstance(node, ast.Name) and node.id == "materials":
            return True
        if isinstance(node, ast.ImportFrom) and any(alias.name == "materials" for alias in node.names):
            return True
    return False


def scanned_modules():
    """The modules (paths relative to `planetgen/`) whose string constants name rows: those that use the table, and `tuning.py`."""
    found = {}
    for path in sorted(SOURCE.rglob("*.py")):
        relative = path.relative_to(SOURCE)
        if relative == TABLE_MODULE:
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"))
        if relative in ALWAYS_SCANNED or _uses_table(tree):
            found[relative] = tree
    return found


def names_of(row):
    return {row.id.lower(), row.name.lower(), *(alias.lower() for alias in row.aliases)}


def reached_rows():
    """Row id -> the first module that names it, for every row a draw reads (a group's members through the group)."""
    named = {}
    for relative, tree in scanned_modules().items():
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and isinstance(node.value, str):
                named.setdefault(node.value.lower(), str(relative))
    reached = {}
    for row in materials.MATERIALS:
        module = next((named[name] for name in sorted(names_of(row)) if name in named), None)
        if module is not None:
            reached[row.id] = module
    for group_id in [row_id for row_id in reached if materials.material(row_id).is_group]:
        for member in materials.expand(group_id):
            reached.setdefault(member.id, f"{reached[group_id]} (through the group {group_id})")
    return reached


def doubles():
    """The sets of rows that are one material under two names: the same formula where a row is a powder or grain-size variant of the rest."""
    by_formula = {}
    for row in materials.MATERIALS:
        if row.formula and (VARIANT_SUFFIX.search(row.id) or not row.crystal_system):   # a row with its own crystal system is a polymorph, not a copy
            by_formula.setdefault(row.formula, []).append(row.id)
    return [frozenset(ids) for ids in by_formula.values() if len(ids) > 1 and any(VARIANT_SUFFIX.search(row_id) for row_id in ids)]


def test_the_reach_scan_looks_at_the_modules_that_use_the_table():
    modules = {str(path) for path in scanned_modules()}
    assert {"physics/ground_albedo.py", "physics/small_body_composition.py", "physics/surface_composition.py", "tuning.py"} <= modules
    reached = reached_rows()
    assert reached["forsterite"] and reached["water_ice"]
    assert "basalt_powder" in reached


def test_a_group_reaches_its_members():
    reached = reached_rows()
    groups = [row for row in materials.MATERIALS if row.is_group and row.id in reached]
    assert groups
    for group in groups:
        assert {member.id for member in materials.expand(group.id)} <= set(reached), group.id


def test_every_row_is_drawn_or_tagged_vocabulary_only():
    reached = reached_rows()
    orphans = [row.id for row in materials.MATERIALS if row.id not in reached and row.id not in materials.VOCABULARY_ONLY]
    assert not orphans, (f"{len(orphans)} rows of the materials table are read by no draw and not tagged: {', '.join(orphans)}. Wire the row to a recipe, or list it in "
                         "VOCABULARY_ONLY in scripts/materials_data.py with why (and run scripts/build_materials_table.py)")


def test_a_tagged_row_is_not_reached():
    reached = reached_rows()
    stale = {row_id: reached[row_id] for row_id in materials.VOCABULARY_ONLY if row_id in reached}
    assert not stale, f"rows tagged vocabulary-only that a draw reads now (remove the tag): {stale}"


def test_every_tag_gives_a_reason_for_a_row_of_the_table():
    ids = {row.id for row in materials.MATERIALS}
    for row_id, reason in materials.VOCABULARY_ONLY.items():
        assert row_id in ids and len(reason.split()) >= 4, row_id


def test_the_doubles_are_the_acknowledged_ones():
    found = set(doubles())
    assert found == set(ACKNOWLEDGED_DOUBLES), (f"doubles found: {sorted(sorted(group) for group in found)}; acknowledged: {sorted(sorted(group) for group in ACKNOWLEDGED_DOUBLES)}. "
                                                "Merge a new double away or name it, with why, in ACKNOWLEDGED_DOUBLES")
    for group, why in ACKNOWLEDGED_DOUBLES.items():
        assert len(group) >= 2 and why


def test_an_orphan_with_no_tag_would_fail(monkeypatch):
    """The gate has teeth: with a tag taken away, the test above names the row."""
    victim = next(iter(materials.VOCABULARY_ONLY))
    trimmed = {key: value for key, value in materials.VOCABULARY_ONLY.items() if key != victim}
    monkeypatch.setattr(materials, "VOCABULARY_ONLY", trimmed)
    try:
        test_every_row_is_drawn_or_tagged_vocabulary_only()
    except AssertionError as error:
        assert victim in str(error)
    else:
        raise AssertionError("an untagged orphan did not fail")
