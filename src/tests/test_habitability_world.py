"""
The habitability score of every planet and moon (GEN.89,
`planetgen.physics.habitability_world`): the world each kind of body maps
to, what generation stores and what the search and the system page show.
"""

import types

import pytest

from planetgen.db import store
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.compact_remnant import NeutronStar
from planetgen.generation.system import StarSystem
from planetgen import tuning
from planetgen.physics import atmosphere, body_radiation, constants, equipment
from planetgen.physics import habitability_world as hw
from planetgen.physics import habitability, life_scores
from planetgen.physics.protection_limits import limit
from planetgen.util import draw

from tests.air_support import fold_air
from tests.fuzz_support import deterministic_entropy
from tests.test_db_persistence import _placed_nebula, _system_point_pc
from tests.test_web_pages import db_client  # noqa: F401 -- fixture
from tests.test_web_search import _panel


def _earth(**changes):
    values = {
        "body_type": "t", "planet_class": "M", "atmospheric_pressure_pa": 101325.0,
        "surface_temperature_k": 288.0, "gravity_g": 1.0, "radius_km": 6371.0, "mass_kg": 5.972e24,
        "p_o2_kpa": 21.2, "p_n2_kpa": 78.0, "p_ar_kpa": 0.9, "p_co2_kpa": 0.04, "p_h2o_kpa": 1.2,
        "hydrosphere": "surface ocean", "ocean_class": "neutral", "ocean_ph": 8.1, "water_activity": 0.98,
        "phosphorus": "limited", "water_mass_fraction": 2.34e-4, "surface_dose_msv_yr": 2.4,
        "dose_ground_msv_yr": 0.5, "energy_flux_w_m2": 82.0, "native_life": "biosphere",
        "life_class": 1, "life_class_variant": "a", "star_teff_k": 5772.0,
    }
    values.update(changes)
    return fold_air(values)


def test_earth_is_blue_everywhere_and_needs_no_kit():
    scored = hw.compute(_earth())
    assert all(scored[name] in (0, None) for name in equipment.FIELDS)  # nothing is needed, so nothing is listed
    assert equipment.line(scored) == "" and equipment.code(scored) == ""
    assert [scored[f"tier_{domain}"] for domain in hw.DOMAINS] == [0] * 6
    assert scored["phi4"] == pytest.approx(1.0)
    assert scored["phi_bio"] > 0.8 and scored["phi_cpx"] > 0.8


def test_the_air_decides_the_kit():
    carbon_dioxide = hw.compute(_earth(p_co2_kpa=3.0))
    assert carbon_dioxide["eq_co2_scrub"] == 2 and carbon_dioxide["eq_supply"] == 0  # a scrubber on the intake
    assert carbon_dioxide["tier_air"] > 0
    thin = hw.compute(_earth(atmospheric_pressure_pa=30000.0, p_o2_kpa=6.0, p_n2_kpa=23.5))
    assert thin["eq_supply"] == 1 and thin["eq_body"] == 0  # a gas-mix mask: enough pressure, too little oxygen
    vacuum = hw.compute(_earth(atmospheric_pressure_pa=0.0, **{f"p_{g}_kpa": 0.0 for g in ("o2", "n2", "ar", "co2", "h2o")}))
    assert vacuum["eq_body"] == 4 and vacuum["eq_supply"] >= 2  # a pressure suit on its own air
    assert vacuum["tier_pressure"] == 3


def test_every_toxic_gas_of_the_table_reaches_the_world():
    """GEN.657: the world reads every gas of the toxic table, not only CO, H2S and SO2 (an ammonia air is a hazard to the filter)."""
    from planetgen.physics import protection_limits
    values = _earth()
    air = dict(values[atmosphere.GAS_FIELD])
    for formula in protection_limits.TOXIC:
        air[formula] = 0.5 * 1000.0
    values[atmosphere.GAS_FIELD] = air
    world = hw.world_from(values)
    for formula in protection_limits.TOXIC:
        assert world.partial_kpa(formula) == pytest.approx(0.5), formula
    clean = hw.compute(_earth())
    foul = hw.compute(values)
    assert foul["eq_gas_filter"] > clean["eq_gas_filter"] or foul["eq_body"] > clean["eq_body"]


def test_a_lethal_dose_needs_a_hardened_habitat():
    scored = hw.compute(_earth(surface_dose_msv_yr=50000.0))
    assert scored["eq_shelter"] == 2
    assert scored["tier_radiation"] == 3


def test_an_ocean_under_ice_sees_only_the_ground_but_has_little_light():
    world = hw.world_from(_earth(hydrosphere="ice-covered ocean", surface_dose_msv_yr=500000.0,
                                 dose_ground_msv_yr=0.5))
    assert world.solvent == "water"
    assert world.habitat_dose_sv_yr == pytest.approx(0.0005)
    assert world.surface_dose_sv_yr == pytest.approx(500.0)
    assert world.water_source == "ice"
    flux = hw.energy_flux_w_m2(_earth(hydrosphere="ice-covered ocean"), 3.8e26, 5.2, 0.05)
    assert flux == pytest.approx(5e-4)


def test_light_is_six_percent_of_the_stars_unless_the_heat_is_more():
    sun = 3.828e26
    flux = hw.energy_flux_w_m2(_earth(), sun, 1.0, 0.09)
    assert flux == pytest.approx(0.06 * 1361, rel=0.02)
    assert hw.energy_flux_w_m2(_earth(), None, None, 0.09) == pytest.approx(9e-4)
    assert hw.energy_flux_w_m2(_earth(), sun, 1000.0, 0.09) == pytest.approx(9e-4)


def test_a_methane_sea_is_a_solvent_but_a_dry_world_is_not():
    titan = _earth(hydrosphere="ice", surface_temperature_k=94.0, atmospheric_pressure_pa=146700.0,
                   p_ch4_kpa=20.0, p_n2_kpa=125.0)
    from planetgen.physics import atmosphere
    saturated = atmosphere.vapour_pressure_kpa("ch4", 94.0)
    titan["gas_partials_pa"]["CH4"] = saturated * 1000.0
    assert hw.solvent(titan) == "hydrocarbon"
    assert hw.world_from(titan).phosphate_umol_l == hw.NON_AQUEOUS_PHOSPHATE_UMOL_L
    dry = _earth(hydrosphere="dry", surface_temperature_k=300.0, water_mass_fraction=1e-7, p_h2o_kpa=0.0)
    assert hw.solvent(dry) is None
    assert hw.compute(dry)["phi_bio"] == 0.0
    assert hw.water_source(dry) is None
    assert hw.water_source(_earth(hydrosphere="dry", water_mass_fraction=1e-4)) == "hydrated"


def test_carbon_and_nitrogen_follow_the_water_or_the_air():
    inventories = hw.inventories(_earth())
    assert inventories["H"] == pytest.approx(1.0, rel=0.05)
    assert inventories["C"] == pytest.approx(1.0, rel=0.05)
    assert inventories["N"] == pytest.approx(1.0, rel=0.1)
    dry_venus = hw.inventories(_earth(water_mass_fraction=1e-7, atmospheric_pressure_pa=9.2e6, p_co2_kpa=8900.0,
                                      p_n2_kpa=300.0, gravity_g=0.9))
    assert dry_venus["H"] < 0.01
    assert dry_venus["C"] > 1.0 and dry_venus["N"] > 1.0


def test_a_chloride_brine_is_chaotropic_and_a_neutral_ocean_is_not():
    brine = hw.world_from(_earth(ocean_class="chloride brine", water_activity=0.5, ocean_ph=5.0))
    assert brine.chaotropicity_kj_kg == pytest.approx(100.0)
    assert hw.world_from(_earth()).chaotropicity_kj_kg == 0.0


def test_humidity_is_the_vapour_against_saturation():
    assert hw.relative_humidity(_earth()) == pytest.approx(1.2 / hw.saturation_kpa(14.85), rel=0.01)
    assert hw.saturation_kpa(0.0) == pytest.approx(0.611, rel=0.01)
    assert 0.05 <= hw.relative_humidity(_earth(p_h2o_kpa=0.0)) <= 0.99


def test_a_gas_giant_scores_nothing():
    assert set(hw.compute({"body_type": "g", "surface_temperature_k": 120.0}).values()) == {None}


def _bodies(system):
    for planet in system.planets:
        if planet.body_type == "a":
            continue
        yield planet
        yield from planet.moons


def test_every_generated_body_is_scored_consistently():
    seen = {"rocky": 0, "giant": 0}
    for seed in range(60):
        if seed >= 12 and seen["giant"]:      # twelve systems, and more only until a gas giant has been seen (whether one is drawn moves with the epochs)
            break
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(SystemConfig(), galactic_center_dist_ly=26000.0)
        for body in _bodies(system):
            if body.body_type == "g":
                seen["giant"] += 1
                assert all(getattr(body, name) is None for name in hw.BODY_FIELDS)
                continue
            seen["rocky"] += 1
            assert 0 <= body.eq_body <= 6 and 0 <= body.eq_supply <= 4
            assert body.energy_flux_w_m2 >= 0.0
            for domain in hw.DOMAINS:
                score, tier = getattr(body, f"phi4_{domain}"), getattr(body, f"tier_{domain}")
                assert 0.0 <= score <= 1.0 and tier in (0, 1, 2, 3)
                assert (tier == 0) == (score >= 1.0)
            for name in ("phi4", "phi_bio", "phi_cpx", "l_solv", "l_chem", "l_ener", "l_rad"):
                assert 0.0 <= getattr(body, name) <= 1.0
            for name in life_scores.FIELDS:
                assert 0 <= getattr(body, name) <= tuning.LIFE_SCORE_MAX
            assert body.phi_cpx <= body.phi_bio + 1e-9
            if body.tier_radiation == 3:
                world = hw.world_from(hw.values_from_body(body))
                if world.dose_sv_yr > habitability.HUMAN.dose[6]:
                    assert body.eq_rad_suit or body.eq_shelter  # a Red dose needs radiation gear
                else:
                    # Red from the ultraviolet alone (GEN.445): the body ladder covers it, not the dose gear
                    assert world.uv_index > limit("uv_cover_index")
                    assert body.eq_body >= (2 if world.uv_index > limit("uv_suit_index") else 1)
    assert seen["rocky"] > 10 and seen["giant"] > 0


def test_a_neutron_stars_planets_get_the_heaviest_radiation_gear():
    seen = 0
    for seed in range(6):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            cfg = SystemConfig()
            system = StarSystem(system_config=cfg, compact_remnant=NeutronStar(cfg))
        for body in _bodies(system):
            if body.body_type == "g":
                continue
            seen += 1
            assert body.eq_shelter == 2 and body.eq_rad_suit == 3 and body.tier_radiation == 3
            assert hw.compact_host_name(body.star) in ("pulsar", "neutron star")
    assert seen


def test_hosts_are_named_by_their_class_and_a_black_hole_planet_is_lethal():
    assert hw.compact_host_name(types.SimpleNamespace(yerkes_class="BH")) == "black hole"
    assert hw.compact_host_name(types.SimpleNamespace(yerkes_class="NS", pulsar_type="young")) == "pulsar"
    assert hw.compact_host_name(types.SimpleNamespace(yerkes_class="NS", pulsar_type="non-pulsing")) == "neutron star"
    assert hw.compact_host_name(types.SimpleNamespace(yerkes_class="V")) is None
    lethal = hw.compute(_earth(), compact_host="black hole")
    assert lethal["eq_shelter"] == 2 and lethal["tier_radiation"] == 3
    assert lethal["phi_bio"] == 0.0


def _system_with_rocky_worlds():
    for seed in range(40):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.BINARY_SYSTEM = False
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(system_config=cfg, galactic_center_dist_ly=26000.0)
        if any(body.body_type == "t" for body in _bodies(system)):
            return cfg, system
    raise AssertionError("no seed gave a rocky world")


def test_the_score_round_trips_and_a_cloud_rescores_the_dose(mysql_config):
    cfg, system = _system_with_rocky_worlds()
    sector = SpaceSector("Scored", edge_ly=11.5)
    sector.add_system(system, position=(0.0, 0.0, 0.0), system_config=cfg)
    sector_id = store.save_sector(sector, config=mysql_config, galaxy_position={
        "center_x_pc": 100.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 100.0,
    })
    conn = store.get_connection(mysql_config)
    try:
        saved = {row["name"]: row for row in conn.execute("SELECT * FROM planets WHERE body_type = 't'").fetchall()}
        checked = 0
        for planet in system.planets:
            if planet.body_type == "t" and planet.name in saved:
                checked += 1
                for name in hw.BODY_FIELDS:
                    if getattr(planet, name) is None:
                        assert saved[planet.name][name] is None, name
                    elif isinstance(getattr(planet, name), str):
                        assert saved[planet.name][name] == getattr(planet, name), name
                    else:
                        assert saved[planet.name][name] == pytest.approx(getattr(planet, name)), name
        assert checked
        with conn:
            # An airless, worst-case world: the cosmic-ray multiplier changes its dose and so its tier.
            conn.execute("UPDATE planets SET atmospheric_pressure_pa = 0, dose_gcr_msv_yr = 500,"
                         " dose_sep_msv_yr = 1, dose_ground_msv_yr = 0.5, surface_dose_msv_yr = 501.5,"
                         " dose_helio_mult = 1, tier_radiation = NULL, eq_body = NULL"
                         " WHERE body_type = 't'")
            cloud = _placed_nebula(conn, sector_id, _system_point_pc(conn, sector_id), radius_ly=3.0, nebula_class="C")
            conn.execute("UPDATE nebulae SET density_cm3 = 3000, temperature_k = 20 WHERE id = ?", (cloud,))
            home = conn.execute("SELECT center_x_pc FROM nebulae WHERE id = ?", (cloud,)).fetchone()["center_x_pc"]
            conn.execute("UPDATE nebulae SET center_x_pc = center_x_pc + 50 WHERE id = ?", (cloud,))
            store.refresh_containment(conn, [sector_id])
            conn.execute("UPDATE nebulae SET center_x_pc = ? WHERE id = ?", (home, cloud))
            store.refresh_containment(conn, [sector_id])
        stored = conn.execute("SELECT * FROM planets WHERE body_type = 't'").fetchall()
        assert stored
        # The re-score reads a row with its air and its star's temperature (`with_life_inputs`), which the life classes need.
        rows = store.with_life_inputs(conn, stored, False)
        for row in rows:
            assert row["dose_helio_mult"] == pytest.approx(2.5, abs=0.01)
            again = hw.compute({**row, "is_moon": False, "native_life": store._stored_native_life(row)})
            assert row["eq_body"] >= 4  # an airless world needs a pressure suit at least: a hot or a cold surface asks for a heavier body
            for name in hw.SCORE_FIELDS:
                assert row[name] == pytest.approx(again[name]), name
    finally:
        conn.close()


def test_the_search_finds_bodies_by_the_kit_they_need(db_client, mysql_config):  # noqa: F811
    cfg, system = _system_with_rocky_worlds()
    sector = SpaceSector("Kitted", edge_ly=11.5)
    sector.add_system(system, position=(0.0, 0.0, 0.0), system_config=cfg)
    store.save_sector(sector, config=mysql_config)
    planets = [p for p in system.planets if p.body_type != "a" and p.eq_body is not None]
    ticks = ("clothing", "thermal_suit", "sealed_suit", "pressure_suit", "armoured_suit", "vehicle_only")
    levels = {p.eq_body for p in planets if p.eq_body}
    assert levels
    level = sorted(levels)[0]
    tick = ticks[level - 1]
    html = db_client.get("/search").get_data(as_text=True)
    assert "Planet: Equipment a Human Needs" in html
    panel = _panel(db_client.get(f"/search?equipment={tick}").get_data(as_text=True), "planets")
    names = {p.name for p in planets if p.eq_body == level}
    assert names and all(name in panel for name in names)
    sample = next(p for p in planets if p.eq_body == level)
    for name in equipment.required_names({name: getattr(sample, name) for name in equipment.FIELDS}):
        assert "flag-eq-required" in panel, panel[panel.index("<tbody"):][:3000]
        assert name[:1].upper() + name[1:] in panel  # UX.154: one chip per item
    others = {p.name for p in planets if p.eq_body != level}
    assert not any(f">{name}<" in panel for name in others)
    # A stray value filters nothing in and does not crash.
    assert db_client.get("/search?equipment=9").status_code in (200, 302)


def test_the_page_chips_are_solid_for_required_and_outlined_for_recommended_in_their_domains_colour():
    from planetgen.web.lib import habitability_chips
    scored = hw.compute(_earth(p_co2_kpa=3.0, surface_dose_msv_yr=500.0))
    body = {**scored, "equipment": equipment.summary(scored)}
    chip = habitability_chips.equipment_chips(body)
    assert "flag-eq-required" in chip and "CO2 scrubber" in chip and "title=" not in chip
    assert "flag-tier-2" in chip  # the CO2 scrubber comes from the Air domain, which is Yellow here
    assert chip.count("flag-eq-recommended") == len(equipment.recommended_names(scored))
    assert habitability_chips.equipment_chips({"equipment": None}) == ""
    assert habitability_chips.equipment_chips({}) == ""


def test_a_molten_or_deep_frozen_world_scores_without_overflowing():
    for kelvin in (3.0, 40.0, 1500.0, 4000.0, 20000.0):
        scored = hw.compute(_earth(surface_temperature_k=kelvin, hydrosphere="dry", p_h2o_kpa=0.0))
        assert scored["tier_temperature"] == 3


def test_harshness_follows_the_worst_domain_and_defaults_to_the_middle():
    assert hw.harshness(hw.compute(_earth())) == 0.0
    cold = hw.compute(_earth(surface_temperature_k=60.0, hydrosphere="dry", p_h2o_kpa=0.0))
    assert hw.harshness(cold) == pytest.approx(1.0 - cold["phi4_temperature"])
    assert hw.harshness(cold) > 0.6  # a Red domain
    assert hw.harshness({}) == hw.DEFAULT_HARSHNESS == 0.5
    gas_giant = hw.compute(_earth(body_type="g"))
    assert hw.harshness(gas_giant) == 0.5
    assert hw.harshness(types.SimpleNamespace(phi4_pressure=1.0, phi4_temperature=0.8, phi4_air=1.0, phi4_water=1.0,
                                              phi4_radiation=0.9, phi4_suitability=0.0)) == pytest.approx(0.2)


# --- GEN.207: Suitability reads the native life -------------------------------------------

def test_the_native_life_is_the_latest_milestone_the_timeline_reached():
    assert hw.native_life(None, "technological_civilization", True) == "none"
    assert hw.native_life("Chlorophyll", None, False) == "microbial"  # life with no timeline (a moon)
    assert hw.native_life("Chlorophyll", None, True) == "none"  # a timeline that has reached no milestone
    for key, stage in (("abiogenesis", "microbial"), ("photosynthesis", "microbial"), ("complex_cells", "microbial"),
                       ("multicellularity", "macroscopic"), ("technological_civilization", "biosphere")):
        assert hw.native_life("Chlorophyll", key, True) == stage, key


def test_suitability_follows_the_native_life_and_microbes_in_the_ocean_under_ice_are_no_food_to_reach():
    for stage, expected in (("none", 0.2), ("microbial", 0.5), ("macroscopic", 1.0), ("biosphere", 1.0)):
        scored = hw.compute(_earth(native_life=stage))
        assert scored["phi4_suitability"] == pytest.approx(expected)
    buried = hw.compute(_earth(hydrosphere="ice-covered ocean", native_life="microbial"))
    assert buried["phi4_suitability"] == pytest.approx(0.2) and buried["tier_suitability"] == 3
    # Boss, 2026-10-24: flora, fauna and a civilization on a world are food even when its ocean is under ice (GEN.765)
    # (the animals of a sub-ice world are a chemical class's, Class 4b here, with Fauna and no Flora: GEN.889)
    farmed = hw.compute(_earth(hydrosphere="ice-covered ocean", native_life="biosphere", life_stage="technological_civilization",
                               evolutionary_pace="Earth Norm", life_class=4, life_class_variant="b", life_complexity=1.0,
                               life_animal_grade=3, internal_heat_flux_w_m2=0.09, surface_temperature_k=268.0))
    assert farmed["life_flora"] == 0 and farmed["life_fauna"] > 0 and farmed["phi4_suitability"] > 0.2


def test_flora_is_always_good_food_and_flora_with_fauna_is_better():
    """Boss, 2026-10-23: the food rung is made from the Flora and Fauna scores."""
    mats = hw.compute(_earth(native_life="microbial", life_stage="photosynthesis", evolutionary_pace="Earth Norm"))
    assert (mats["life_flora"] > 0, mats["life_fauna"]) == (True, 0)
    assert mats["phi4_suitability"] == pytest.approx(0.85)  # flora alone: Green
    microbes = hw.compute(_earth(native_life="microbial", life_stage="abiogenesis", evolutionary_pace="Earth Norm"))
    assert (microbes["life_microbial"] > 0, microbes["life_flora"]) == (True, 0)
    assert microbes["phi4_suitability"] == pytest.approx(0.5)
    both = hw.compute(_earth(native_life="macroscopic", life_stage="multicellularity", evolutionary_pace="Earth Norm"))
    assert both["life_flora"] > 0 and both["life_fauna"] > 0
    assert both["phi4_suitability"] == pytest.approx(1.0)
    assert life_scores.food_stage(0, 0, 3) == "macroscopic"  # fauna is all animal life: food on its own


def test_a_dry_world_has_red_water_and_red_suitability_and_an_overall_at_most_yellow():
    dry = hw.compute(_earth(hydrosphere="dry", ocean_ph=None, water_activity=None, p_h2o_kpa=0.0,
                            water_mass_fraction=0.0, native_life="biosphere"))
    assert dry["tier_water"] == 3 and dry["phi4_water"] == 0.0
    assert dry["tier_suitability"] == 3
    assert hw.TIER_INDEX["Yellow"] == 2 and 0.4 <= dry["phi4"] < 0.75


def test_life_is_drawn_after_the_score_and_the_score_is_redone_with_it():
    """Suitability never beats the food rung of the life the body was given, and some bodies have life."""
    stages = set()
    for seed in range(12):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            config = SystemConfig()
            config.STAR_TYPE = "G2V"  # a cool star's zone has lost its water before any life (GEN.639)
            system = StarSystem(config, galactic_center_dist_ly=26000.0)
        for body in _bodies(system):
            if body.body_type == "g":
                continue
            stage = life_scores.food_stage(body.life_microbial, body.life_flora, body.life_fauna)
            stages.add(stage)
            food = hw.habitability.FOOD_SCORE["none" if body.hydrosphere == "ice-covered ocean" and stage == "microbial" else stage]
            assert body.phi4_suitability <= food + 1e-9, (body.name, stage)
            assert hw.TIER_INDEX["Red"] == body.tier_suitability or body.phi4_suitability >= 0.4
    assert stages - {"none"}, stages


def test_the_equipment_cell_draws_one_chip_per_item_in_its_domains_colour():
    """UX.154: the search cell carries `chips` (solid when required, outlined when recommended) and a muted dash for none."""
    from planetgen.web.lib import habitability_cells
    row = {"equipment": "Required: gas-mix mask. Recommended: hardened habitat.", "tier_air": 2, "tier_radiation": 3,
           "equipment_items": [{"item": "gas-mix mask", "severity": "required", "domain": "Air"},
                               {"item": "hardened habitat", "severity": "recommended", "domain": "Radiation"}]}
    cell = habitability_cells.equipment_cell(row)
    assert cell["chips"] == [{"text": "gas-mix mask", "severity": "required", "tier": 2, "anchor": "equipment-gas-mix-mask"},
                             {"text": "hardened habitat", "severity": "recommended", "tier": 3,
                              "anchor": "equipment-shelter"}]
    assert cell["text"] == row["equipment"]
    assert habitability_cells.equipment_cell({"equipment": "", "equipment_items": []}) == {"text": "—", "muted": True}


def _sun_body(**fields):
    base = {"planet_class": "M", "gravity": 1.0, "atmospheric_pressure": 101325.0, "distance": 1.0, "is_moon": False,
            "star": types.SimpleNamespace(temperature=5772.0, luminosity=constants.SOLAR_LUMINOSITY, uv_c_share=1.0),
            "p_n2_kpa": 79.1, "p_o2_kpa": 21.2, "p_ar_kpa": 0.95, "ozone_column_du": 300.0}
    base.update(fields)
    return types.SimpleNamespace(**fold_air(base))


def test_the_ozone_column_moves_the_radiation_tier_through_the_ground_ultraviolet():
    """GEN.445: the habitability reads the derived UV, so changing the ozone of an otherwise Earth-like world moves the Radiation tier (and
    no other domain)."""
    def scored(ozone_du):
        index = body_radiation.dna_dose_index(_sun_body(ozone_column_du=ozone_du))
        return hw.compute(_earth(uv_surface_index=index))

    shielded, thin, bare = scored(300.0), scored(60.0), scored(0.0)
    assert shielded["tier_radiation"] == 0
    assert bare["tier_radiation"] > shielded["tier_radiation"] and bare["tier_radiation"] >= thin["tier_radiation"]
    assert bare["tier_radiation"] == 3                    # the unfiltered light of a world with oxygen and no ozone is lethal
    for domain in ("pressure", "temperature", "air", "water"):
        assert bare[f"tier_{domain}"] == shielded[f"tier_{domain}"]


def test_the_habitability_cell_carries_the_four_life_scores_after_the_phi6_symbols():
    """UX.152: a scored body's cell has its life cluster, a body without scores has none."""
    from planetgen.web.lib import habitability_cells
    row = {"phi4": 0.9, "tier_pressure": 0, "tier_temperature": 0, "tier_air": 0, "tier_water": 0, "tier_radiation": 0,
           "tier_suitability": 0, "life_microbial": 3, "life_flora": 0, "life_fauna": 1, "life_sentient": 0}
    cell = habitability_cells.habitability_cell(row)
    assert [entry["name"] for entry in cell["life"]] == ["microbial", "flora", "fauna", "sentient"]
    assert [entry["score"] for entry in cell["life"]] == [3, 0, 1, 0]
    assert cell["life"][1]["label"] == "Flora: none" and cell["life"][0]["anchor"] == "life-microbial"
    assert "Microbial life: 3 of 9" in cell["text"]
    assert "life" not in habitability_cells.habitability_cell({**row, "life_flora": None})
