# tests/test_phenomena_registry.py

"""
GEN.250: the phenomena registry (`planetgen/generation/phenomena/registry.py`)
is the one table of every kind of object. These tests check that every row
is complete, that everything the registry names resolves, and that every
other place in the code that names a kind agrees with it, so a kind cannot
exist without a rate, a builder, an ID code and a page.
"""

import ast
import pathlib
import re

import pytest

from planetgen import tuning
from planetgen.db import models, store
from planetgen.galaxy import sector as galaxy_sector
from planetgen.generation import phenomenon_scatter, run_phenomenon, run_sector
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena import hosts, registry
from planetgen.names import naming_key, object_id
from planetgen.web.lib import classref

_SOURCE = pathlib.Path(registry.__file__).resolve().parents[2]
_NOT_PHENOMENON_PAGES = {"star-spectral", "star-luminosity", "planet"}


def test_every_row_has_every_column_and_valid_values():
    assert tuple(registry.Kind.__dataclass_fields__) == registry.COLUMNS
    assert len({row.kind for row in registry.KINDS}) == len(registry.KINDS)
    for row in registry.KINDS:
        assert row.kind and row.evidence_note and row.status == "built", row.kind
        assert row.family in registry.FAMILIES, row.kind
        assert row.gate in registry.GATES, row.kind
        assert row.mode in registry.MODES, row.kind
        assert row.host is None or row.host in hosts.BY_NAME, row.kind
        assert isinstance(row.carries, tuple) and isinstance(row.sector_roll, bool), row.kind
        assert (row.scatter_pass in registry.SCATTER_PASSES) == (row.mode == "scatter"), row.kind


def test_every_builder_resolves_and_a_made_kind_has_one():
    for row in registry.KINDS:
        if row.mode in ("scatter", "fill", "field", "flag") or row.sector_roll:
            assert row.builder, f"{row.kind} has no builder"
        if row.builder:
            assert callable(registry.resolve(row.builder)), f"{row.kind}: {row.builder}"


def test_every_sector_roll_builder_takes_a_config_and_a_distance():
    for kind, _type, factory in run_sector.SECTOR_PHENOMENON_KINDS:
        assert factory(SystemConfig(), 1000.0) is not None, kind


def test_carried_kinds_and_page_classes_resolve():
    for row in registry.KINDS:
        for carried in row.carries:
            assert carried in registry.BY_KIND, f"{row.kind} carries {carried}"
    slugs = set(classref.catalog()) - _NOT_PHENOMENON_PAGES
    pages = {row.page_class for row in registry.KINDS if row.page_class}
    assert pages == slugs


def test_no_id_code_repeats_and_codes_fit_six_bits():
    codes = [row.id_code for row in registry.KINDS if row.id_code is not None]
    assert len(codes) == len(set(codes)) and all(1 <= code < 64 for code in codes)


def test_the_existing_lists_are_the_registry():
    assert object_id.KIND_CODES == registry.id_codes()
    assert set(naming_key.KINDS) == set(registry.id_codes()) - {"bright-star"} | {"constellation"}
    assert phenomenon_scatter.SCATTERED_KINDS == registry.kinds(scatter_pass="band-density")
    assert tuple(kind for kind, _t, _f in run_sector.SECTOR_PHENOMENON_KINDS) == registry.kinds(sector_roll=True)
    assert tuning.PHENOMENON_TYPE_CHOICES == registry.phenomenon_types()


def test_the_density_table_and_the_registry_name_the_same_kinds():
    assert set(tuning.PHENOMENON_DENSITY_PC3) == set(registry.rated_kinds())
    assert set(tuning.PHENOMENON_RATE_SCALE) == set(tuning.PHENOMENON_DENSITY_PC3)
    for row in registry.KINDS:
        assert row.rate in (None, row.kind), row.kind
    # A kind the sector rolls has a rate (or is explicitly rate-less: the black holes' follows the stars' metallicity,
    # `galaxy.remnant_density`, GEN.563).
    assert {row.kind for row in registry.rows(sector_roll=True) if not row.rate} == {"black-hole"}


def test_no_place_names_a_kind_the_registry_lacks():
    types = set(registry.phenomenon_types())
    id_kinds = set(registry.id_codes())
    assert set(galaxy_sector._PHENOMENON_CLASSES_BY_TYPE) == types
    assert set(store._PHENOMENON_INSERTERS) == types
    assert set(run_phenomenon.TYPE_LABELS) == types
    assert set(store.OBJECT_ID_TABLES) == id_kinds - {"bright-star"}
    assert {kind for kind, _subtype in phenomenon_scatter.SCATTER_CLASSES} <= set(phenomenon_scatter.SCATTERED_KINDS)
    for row in registry.KINDS:
        assert row.id_code is not None or row.phenomenon_type is not None or row.mode != "scatter" \
            or row.scatter_pass in ("bright-star", "wanderer"), row.kind


def test_the_scatter_table_allows_exactly_the_scattered_kinds():
    table = models.metadata.tables["phenomenon_scatter"]
    check = next(c for c in table.constraints if getattr(c, "name", None) == "chk_phenomenon_scatter_kind")
    allowed = set(re.findall(r"'([^']+)'", str(check.sqltext)))
    placed = {row.kind for row in registry.rows(mode="scatter") if row.scatter_pass != "bright-star"}
    assert allowed == placed


def _string_arguments(call, index):
    return [call.args[index].value] if len(call.args) > index and isinstance(call.args[index], ast.Constant) \
        and isinstance(call.args[index].value, str) else []


def test_every_literal_kind_passed_to_add_phenomenon_is_in_the_registry():
    """No code path makes an object of a kind that is not in the table (the
    kind is the second argument of `add_phenomenon`)."""
    types = set(registry.phenomenon_types())
    seen = []
    for path in _SOURCE.rglob("*.py"):
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, ast.Call) and getattr(node.func, "attr", None) == "add_phenomenon":
                seen += _string_arguments(node, 1)
    assert seen, "the scan found no add_phenomenon call"
    assert set(seen) <= types
    for kind in set(seen):
        assert kind in registry.BY_KIND


@pytest.mark.parametrize("kind", [row.kind for row in registry.rows(sector_roll=True)])
def test_a_sector_rolled_kind_is_saved_as_its_phenomenon_type(kind):
    row = registry.BY_KIND[kind]
    assert row.phenomenon_type in galaxy_sector._PHENOMENON_CLASSES_BY_TYPE
