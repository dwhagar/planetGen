"""GEN.427: the shared materials and components table (`planetgen/physics/materials.py`, its data file
`physics/data/materials.json`, built by `scripts/build_materials_table.py` from part E's table, parts B and C and the values in
`scripts/materials_data.py`). The audit fails when a material named in a description or in the class tables is not a row, and
when a number has no source."""

import json
import os

import pytest

from planetgen import tuning
from planetgen.generation import possibility_table as pt
from planetgen.physics import elements, gases, materials
from planetgen.physics.elements import DataError

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")


# --- The file states where each number comes from ---------------------------------------------------------------

def test_every_row_has_a_source_a_state_and_a_gate():
    assert len(materials.MATERIALS) >= 150
    for row in materials.MATERIALS:
        assert row.source and row.state in pt.STATES, row.id
        assert row.gate and (row.gate == (pt.ANY_CLASS,) or set(row.gate) <= set(tuning.PLANET_CLASSES)), row.id
        if row.state != pt.VERIFIED:
            assert row.limit, f"{row.id} is {row.state} and names no limit"


def test_every_number_has_a_source_and_a_code():
    for row in materials.MATERIALS:
        for key, prop in row.props.items():
            if prop.value is None:
                assert prop.code == "N" and prop.source is None, (row.id, key)
            else:
                assert prop.source in materials.SOURCES and prop.code in ("R", "C", "H", "S"), (row.id, key)
        for entry in row.eos:
            assert entry.source in materials.SOURCES and entry.code in ("R", "C", "H", "S"), (row.id, entry.model)


def test_no_recalled_value_is_kept():
    """Part E's RECALLED (Q) cells are dropped to null, not kept: sulfur's melting point, the rock rows' solidus, latent heats."""
    assert materials.material("iron").prop("t_melt_k").value is None
    assert "recalled" in materials.material("iron").note
    assert materials.material("basalt_powder").prop("t_melt_k").value is None
    for row in materials.MATERIALS:
        assert not any(prop.code == "Q" for prop in row.props.values()), row.id


def test_the_labels_name_a_column_of_the_row_with_the_same_code():
    for row in materials.MATERIALS:
        for pair in filter(None, row.labels.split(";")):
            key, _, code = pair.partition("=")
            columns = [k for k, (_unit, label_key) in materials.PROPERTIES.items() if label_key == key]
            assert columns, (row.id, pair)
            assert any(k in row.props and row.props[k].code == code for k in columns), (row.id, pair)


def test_the_ids_and_names_are_unique_and_a_name_finds_its_row():
    ids = [row.id for row in materials.MATERIALS]
    assert len(ids) == len(set(ids))
    for row in materials.MATERIALS:
        assert materials.material(row.id) is row and materials.material(row.name) is row
        for alias in row.aliases:
            assert materials.material(alias.upper()) is row


def test_the_molar_masses_come_from_the_element_table():
    assert materials.material("forsterite").molar_mass == pytest.approx(elements.molar_mass("Mg2SiO4"), rel=0, abs=0)
    assert materials.material("gypsum").molar_mass == pytest.approx(0.17217, abs=2e-5)
    assert materials.material("kamacite").molar_mass is None
    for row in materials.MATERIALS:
        assert (row.molar_mass is None) == (row.formula is None), row.id


def test_part_e_values_are_kept_with_their_labels_and_the_listed_albedos():
    ice = materials.material("water_ice")
    assert ice.prop("p_crit_bar").value == 220.64 and ice.prop("alb_040_050").value == 0.832
    assert ice.prop("rho_grain_g_cm3").value == pytest.approx(0.917)
    assert "Tm=H" in ice.labels
    assert materials.material("komatiite").props == {}


# --- The equations of state --------------------------------------------------------------------------------------

def test_seager_2007_table_1_is_in_the_rows():
    """Read from the paper: epsilon iron 156.2 GPa and 6.08 Vinet, MgO 177.0, MgSiO3 perovskite 247 with the fourth-order K''."""
    iron = next(e for e in materials.material("iron_epsilon").eos if e.source == "seager2007" and e.model == "vinet")
    assert (iron.params["rho0_kg_m3"], iron.params["k0_gpa"], iron.params["k0p"]) == (8300.0, 156.2, 6.08)
    assert iron.uncertainty == {"k0_gpa": 1.8, "k0p": 0.12}
    pv = next(e for e in materials.material("mg_bridgmanite").eos if e.source == "seager2007" and e.model == "bme")
    assert pv.params["order"] == 4 and pv.params["k0pp_per_gpa"] == -0.016 and pv.params["k0_gpa"] == 247.0
    mgo = next(e for e in materials.material("periclase").eos if e.source == "seager2007")
    assert mgo.params["k0_gpa"] == 177.0
    slb = materials.material("periclase").eos_of("mgd-slb3")
    assert slb.params["k0_gpa"] == 161.14, "both periclase choices are kept, each with its source"
    ice_vii = materials.material("ice_vii").eos_of("bme")
    assert (ice_vii.params["k0_gpa"], ice_vii.params["k0p"], ice_vii.params["rho0_kg_m3"]) == (23.7, 4.15, 1460.0)


def test_seager_table_2_is_monotone_and_ends_at_7686_gpa():
    table = materials.material("ice_x").eos_of("tabulated")
    pressures = [p for p, _rho in table.params["points"]]
    densities = [rho for _p, rho in table.params["points"]]
    assert pressures == sorted(pressures) and densities == sorted(densities)
    assert len(pressures) == 26 and pressures[-1] == pytest.approx(7686.171) and densities[0] == pytest.approx(1636.617)


def test_the_polytrope_fits_reproduce_their_zero_pressure_density_and_grow():
    for row_id in ("iron_epsilon", "mg_bridgmanite", "ice_x", "moissanite", "carbon_graphite"):
        fit = materials.material(row_id).eos_of("modified-polytrope")
        rho0, c, n = fit.params["rho0_kg_m3"], fit.params["c_kg_m3_pa_n"], fit.params["n"]
        assert rho0 + c * 1e9 ** n > rho0 and fit.valid["p_pa"][1] == 1e16


def test_every_slb_unit_has_the_molar_mass_of_its_formula():
    count = 0
    for row in materials.MATERIALS:
        for entry in row.eos:
            if entry.model == "mgd-slb3":
                count += 1
                assert abs(elements.molar_mass(entry.params["unit_formula"]) * 1000.0 - entry.params["unit_molar_mass_g_mol"]) < 0.05
    assert count == 34


def test_the_ice_phases_cover_the_pressure_sequence():
    for row_id, low in (("water_ice", 0.1), ("ice_iii", 209.5), ("ice_v", 355.0), ("ice_vi", 618.4), ("ice_vii", 2216.0),
                        ("ice_x", 50000.0)):
        points = materials.material(row_id).eos_of("reference-points").params["points"]
        assert points[0][0] == low
    assert materials.material("water_liquid").prop("rho_grain_g_cm3").value == pytest.approx(0.99656)
    # the ice phases above 1.7 GPa are planetary only; ice Ih is every class
    assert materials.material("ice_vii").gate == materials.material("ice_x").gate == tuple(sorted(tuning.PLANET_CLASSES))
    assert materials.material("water_ice").gate == (pt.ANY_CLASS,)


def test_hydrogen_and_helium_hold_no_invented_constants():
    for name in ("hydrogen", "helium"):
        row = materials.material(name)
        assert [e.model for e in row.eos] == ["ideal-gas"] and row.gate == tuple(sorted(tuning.PLANET_CLASSES))
        assert row.prop("rho_gas_stp_g_l").value > 0


def test_densities_are_the_values_read_from_the_sources():
    assert materials.material("troilite").prop("rho_grain_g_cm3").value == 4.85
    assert materials.material("troilite").prop("rho_grain_g_cm3").measured_range == (4.67, 4.79)
    assert materials.material("quartz").prop("rho_grain_g_cm3").value == pytest.approx(2.66)
    assert materials.material("co2_ice").prop("rho_grain_g_cm3").source == "pubchem"
    assert materials.material("naphthalene").prop("rho_grain_g_cm3").value == pytest.approx(1.162)
    assert materials.material("ch4_ice").prop("rho_grain_g_cm3").value is None, "no solid methane density was reached"
    for name in ("magnesium", "iridium", "gold"):
        assert materials.material(name).prop("rho_grain_g_cm3").source == "crc-density-wiki"


def test_the_gas_pointers_name_rows_of_the_gas_table():
    pointers = [row for row in materials.MATERIALS if row.gas_formula]
    assert len(pointers) >= 10
    for row in pointers:
        assert gases.has(row.gas_formula), (row.id, row.gas_formula)


# --- The groups -------------------------------------------------------------------------------------------------

def test_a_group_has_members_and_no_number_of_its_own():
    groups = [row for row in materials.MATERIALS if row.is_group]
    assert {g.id for g in groups} >= {"olivine", "pyroxene", "plagioclase", "silicate_dust", "hydrocarbons", "rock"}
    for group in groups:
        assert not group.props and not group.eos and group.formula is None
        assert all(not member.is_group for member in materials.expand(group.id))
    assert [row.id for row in materials.expand("olivine")] == ["forsterite", "fayalite"]
    assert [row.id for row in materials.expand("forsterite")] == ["forsterite"]


# --- The audit: every material the lists name is a row (GEN.427) -------------------------------------------------

def test_every_asteroid_component_is_a_row():
    assert len(tuning.ASTEROID_COMPONENTS) == 84
    missing = [name for name in tuning.ASTEROID_COMPONENTS if not materials.has(name)]
    assert not missing, f"not in materials.json: {missing}"


def test_every_field_family_and_comet_component_is_a_row():
    names = list(tuning.COMET_COMPOSITION)
    for family in tuning.ASTEROID_FIELD_COMPOSITIONS.values():
        names += family["components"] or []  # None: any asteroid component, already covered
    missing = sorted({name for name in names if not materials.has(name)})
    assert not missing, f"not in materials.json: {missing}"


def test_every_class_composition_text_is_mapped_to_rows():
    texts = {data["composition"] for data in tuning.PLANET_CLASSES.values()}
    assert len(texts) == 11
    for text in texts:
        rows = materials.class_composition(text)
        assert rows and all(isinstance(row, materials.Material) for row in rows), text
    with pytest.raises(materials.UnknownMaterial, match="not mapped"):
        materials.class_composition("a made-up composition")


def test_an_unknown_material_is_refused_with_the_way_to_add_it():
    with pytest.raises(materials.UnknownMaterial, match="materials.json"):
        materials.material("unobtainium")
    assert not materials.has("unobtainium")


def test_there_is_no_class_list_apart_from_the_shared_table():
    """The materials are registered as the one shared components table that every other table's fractions refer to."""
    table = pt.REGISTRY.get(pt.SHARED_COMPONENTS)
    assert table is materials.TABLE and len(table) == len(materials.MATERIALS)
    pt.REGISTRY.validate()
    assert table.state_counts()[pt.OPEN_STATE] == 0
    # a fraction naming a material the table lacks is refused by the registry
    registry = pt.Registry()
    registry.register(materials.TABLE.__class__(pt.SHARED_COMPONENTS, 1, list(materials.TABLE.rows)))
    registry.register(pt.Table("belts", 1, [pt.Row("a", (pt.ANY_CLASS,), 1.0, "test", fractions={"unobtainium": 1.0})],
                               components=pt.SHARED_COMPONENTS))
    with pytest.raises(pt.TableError, match="unobtainium"):
        registry.validate()


# --- The loader refuses a file or row without what it must state ------------------------------------------------------

@pytest.fixture
def broken(tmp_path, monkeypatch):
    """Returns `edit(function)`: copies the data file with `function(data)` applied and makes the loader read it."""
    original = materials.DATA_FILE

    def edit(function):
        with open(original, encoding="utf-8") as handle:
            data = json.load(handle)
        function(data)
        path = tmp_path / "materials.json"
        path.write_text(json.dumps(data), encoding="utf-8")
        monkeypatch.setattr(materials, "DATA_FILE", str(path))
        return materials._load
    return edit


def _row(data, row_id):
    return next(item for item in data["materials"] if item["id"] == row_id)


@pytest.mark.parametrize("field", ["sources", "retrieved", "title"])
def test_the_loader_refuses_a_file_that_states_no_source(broken, field):
    with pytest.raises(DataError, match=field):
        broken(lambda data: data.pop(field))()


@pytest.mark.parametrize("field", ["source", "state", "gate"])
def test_the_loader_refuses_a_row_without_a_source_state_or_gate(broken, field):
    with pytest.raises(DataError, match=field):
        broken(lambda data: _row(data, "forsterite").pop(field))()


def test_the_loader_refuses_a_number_with_no_source(broken):
    with pytest.raises(DataError, match="no source"):
        broken(lambda data: _row(data, "forsterite")["props"]["rho_grain_g_cm3"].update(source=None))()
    with pytest.raises(DataError, match="no source"):
        broken(lambda data: _row(data, "forsterite")["props"]["rho_grain_g_cm3"].update(source="a rumour"))()


def test_the_loader_refuses_a_value_with_the_code_not_reached_and_a_null_with_another(broken):
    with pytest.raises(DataError, match="code N"):
        broken(lambda data: _row(data, "forsterite")["props"]["rho_grain_g_cm3"].update(code="N"))()
    with pytest.raises(DataError, match="no value"):
        broken(lambda data: _row(data, "forsterite")["props"]["rho_grain_g_cm3"].update(value=None))()


def test_the_loader_refuses_a_partly_verified_row_without_its_limit(broken):
    def edit(data):
        _row(data, "forsterite")["limit"] = ""
    with pytest.raises(DataError, match="names no limit"):
        broken(edit)()


def test_the_loader_refuses_an_unreadable_formula_a_wrong_unit_and_an_unknown_gate(broken):
    with pytest.raises(DataError, match="formula"):
        broken(lambda data: _row(data, "forsterite").update(formula="Qq9"))()
    with pytest.raises(DataError, match="g/cm3"):
        broken(lambda data: _row(data, "forsterite")["props"]["rho_grain_g_cm3"].update(unit="kg/m3"))()
    with pytest.raises(DataError, match="gate"):
        broken(lambda data: _row(data, "forsterite").update(gate="giants"))()


def test_the_loader_refuses_an_equation_of_state_without_the_constants_of_its_model(broken):
    def edit(data):
        del _row(data, "iron_epsilon")["eos"][0]["params"]["k0p"]
    with pytest.raises(DataError, match="k0p"):
        broken(edit)()
    with pytest.raises(DataError, match="unknown model"):
        broken(lambda data: _row(data, "iron_epsilon")["eos"][0].update(model="magic"))()


def test_the_loader_refuses_a_duplicate_id_an_alias_that_is_another_row_and_a_missing_member(broken):
    with pytest.raises(DataError, match="twice"):
        broken(lambda data: _row(data, "fayalite").update(id="forsterite"))()
    with pytest.raises(DataError, match="both"):
        broken(lambda data: _row(data, "fayalite")["aliases"].append("forsterite"))()
    with pytest.raises(DataError, match="not a row"):
        broken(lambda data: _row(data, "olivine").update(members=["forsterite", "nothing"]))()


def test_the_loader_refuses_a_label_that_names_no_column(broken):
    with pytest.raises(DataError, match="labels"):
        broken(lambda data: _row(data, "forsterite").update(labels="rho=R;nonsense=H"))()


def test_the_loader_refuses_a_group_with_a_number_of_its_own(broken):
    with pytest.raises(DataError, match="group"):
        broken(lambda data: _row(data, "olivine").update(formula="Mg2SiO4"))()


# --- The builder -------------------------------------------------------------------------------------------------

def _builder():
    import importlib.util
    spec = importlib.util.spec_from_file_location("build_materials_table", os.path.join(REPO, "scripts", "build_materials_table.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_licence_is_not_a_column_and_has_its_own_file_in_the_root():
    with open(materials.DATA_FILE, encoding="utf-8") as handle:
        text = handle.read()
    assert '"licence"' not in text
    assert not hasattr(materials.material("forsterite"), "licence")
    with open(os.path.join(REPO, materials.LICENCE_FILE), encoding="utf-8") as handle:
        notes = handle.read()
    for key in materials.SOURCES:
        assert f"`{key}`" in notes, f"{materials.LICENCE_FILE} says nothing about the source {key}"


def test_the_properties_of_the_builder_and_the_loader_agree():
    assert _builder().PROPERTIES == materials.PROPERTIES


def test_the_data_file_is_what_the_builder_makes(tmp_path, monkeypatch, capsys):
    """The checked-in file is the research tables' values and nothing else (`scripts/build_materials_table.py --check`)."""
    module = _builder()
    if not os.path.isfile(module.PART_E):
        pytest.skip("the research parts are not in this checkout")
    monkeypatch.setattr("sys.argv", ["build_materials_table.py", "--check"])
    module.main()
    assert "current" in capsys.readouterr().out
