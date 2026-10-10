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
from planetgen.physics import habitability_world as hw
from planetgen.util import draw

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
        "dose_ground_msv_yr": 0.5, "energy_flux_w_m2": 82.0,
    }
    values.update(changes)
    return values


def test_earth_is_blue_everywhere_and_needs_no_kit():
    scored = hw.compute(_earth())
    assert scored["equipment_tier"] == 0
    assert [scored[f"tier_{domain}"] for domain in hw.DOMAINS] == [0, 0, 0, 0]
    assert scored["phi4"] == pytest.approx(1.0)
    assert scored["phi_bio"] > 0.8 and scored["phi_cpx"] > 0.8 and scored["phi_tech"] > 0.95
    assert scored["hab_note"] == "All four domains Blue"


def test_the_air_decides_the_kit():
    carbon_dioxide = hw.compute(_earth(p_co2_kpa=3.0))
    assert carbon_dioxide["equipment_tier"] == 2  # mask and scrubber
    assert "Chemistry" in carbon_dioxide["hab_note"]
    thin = hw.compute(_earth(atmospheric_pressure_pa=30000.0, p_o2_kpa=6.0, p_n2_kpa=23.5))
    assert thin["equipment_tier"] == 1  # a mask: enough pressure, too little oxygen
    vacuum = hw.compute(_earth(atmospheric_pressure_pa=0.0, **{f"p_{g}_kpa": 0.0 for g in ("o2", "n2", "ar", "co2", "h2o")}))
    assert vacuum["equipment_tier"] >= 3  # a sealed suit at best
    assert vacuum["tier_pressure"] == 3


def test_a_lethal_dose_needs_full_life_support():
    scored = hw.compute(_earth(surface_dose_msv_yr=50000.0))
    assert scored["equipment_tier"] == 4
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
    titan["p_ch4_kpa"] = saturated
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
    for seed in range(12):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(SystemConfig(), galactic_center_dist_ly=26000.0)
        for body in _bodies(system):
            if body.body_type == "g":
                seen["giant"] += 1
                assert all(getattr(body, name) is None for name in hw.BODY_FIELDS)
                continue
            seen["rocky"] += 1
            assert 0 <= body.equipment_tier <= 4
            assert body.energy_flux_w_m2 >= 0.0
            for domain in hw.DOMAINS:
                score, tier = getattr(body, f"phi4_{domain}"), getattr(body, f"tier_{domain}")
                assert 0.0 <= score <= 1.0 and tier in (0, 1, 2, 3)
                assert (tier == 0) == (score >= 1.0)
            for name in ("phi4", "phi_bio", "phi_cpx", "phi_tech", "l_solv", "l_chem", "l_ener", "l_rad"):
                assert 0.0 <= getattr(body, name) <= 1.0
            assert body.phi_cpx <= body.phi_bio + 1e-9
            if body.tier_radiation == 3:
                assert body.equipment_tier == 4
    assert seen["rocky"] > 10 and seen["giant"] > 0


def test_a_neutron_stars_planets_get_the_lowest_tier():
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
            assert body.equipment_tier == 4 and body.tier_radiation == 3
            assert body.hab_note.startswith(("Planet of a pulsar", "Planet of a neutron star"))
    assert seen


def test_hosts_are_named_by_their_class_and_a_black_hole_planet_is_lethal():
    assert hw.compact_host_name(types.SimpleNamespace(yerkes_class="BH")) == "black hole"
    assert hw.compact_host_name(types.SimpleNamespace(yerkes_class="NS", pulsar_type="young")) == "pulsar"
    assert hw.compact_host_name(types.SimpleNamespace(yerkes_class="NS", pulsar_type="non-pulsing")) == "neutron star"
    assert hw.compact_host_name(types.SimpleNamespace(yerkes_class="V")) is None
    lethal = hw.compute(_earth(), compact_host="black hole")
    assert lethal["equipment_tier"] == 4 and lethal["tier_radiation"] == 3
    assert lethal["phi_bio"] == 0.0 and "black hole" in lethal["hab_note"]


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
                         " dose_helio_mult = 1, tier_radiation = NULL, equipment_tier = NULL"
                         " WHERE body_type = 't'")
            cloud = _placed_nebula(conn, sector_id, _system_point_pc(conn, sector_id), radius_ly=3.0, nebula_class="C")
            conn.execute("UPDATE nebulae SET density_cm3 = 3000, temperature_k = 20 WHERE id = ?", (cloud,))
            home = conn.execute("SELECT center_x_pc FROM nebulae WHERE id = ?", (cloud,)).fetchone()["center_x_pc"]
            conn.execute("UPDATE nebulae SET center_x_pc = center_x_pc + 50 WHERE id = ?", (cloud,))
            store.refresh_containment(conn, [sector_id])
            conn.execute("UPDATE nebulae SET center_x_pc = ? WHERE id = ?", (home, cloud))
            store.refresh_containment(conn, [sector_id])
        rows = conn.execute("SELECT * FROM planets WHERE body_type = 't'").fetchall()
        assert rows
        for row in rows:
            assert row["dose_helio_mult"] == pytest.approx(2.5, abs=0.01)
            again = hw.compute(dict(row))
            assert row["equipment_tier"] is not None and row["equipment_tier"] >= 3  # an airless world
            for name in hw.SCORE_FIELDS:
                assert row[name] == pytest.approx(again[name]), name
    finally:
        conn.close()


def test_the_search_finds_bodies_by_the_kit_they_need(db_client, mysql_config):  # noqa: F811
    cfg, system = _system_with_rocky_worlds()
    sector = SpaceSector("Kitted", edge_ly=11.5)
    sector.add_system(system, position=(0.0, 0.0, 0.0), system_config=cfg)
    store.save_sector(sector, config=mysql_config)
    tiers = {planet.equipment_tier for planet in system.planets
             if planet.body_type != "a" and planet.equipment_tier is not None}
    assert tiers
    tier = sorted(tiers)[0]
    html = db_client.get("/search").get_data(as_text=True)
    assert "Planet: Equipment a Human Needs" in html
    panel = _panel(db_client.get(f"/search?equipment={tier}").get_data(as_text=True), "planets")
    planets = [p for p in system.planets if p.body_type != "a"]
    names = {p.name for p in planets if p.equipment_tier == tier}
    assert names and all(name in panel for name in names)
    assert hw.EQUIPMENT_LABELS[tier] in panel
    others = {p.name for p in planets if p.equipment_tier is not None and p.equipment_tier != tier}
    assert not any(f">{name}<" in panel for name in others)
    # A stray value filters nothing in and does not crash.
    assert db_client.get("/search?equipment=9").status_code in (200, 302)


def test_the_page_chip_names_the_kit_and_takes_the_worst_domains_colour():
    from planetgen.web.lib import systempage
    scored = hw.compute(_earth(p_co2_kpa=3.0, surface_dose_msv_yr=5000.0))
    chip = systempage._habitability_chip(scored)
    assert "flag-tier-2" in chip and "Mask and scrubber" in chip and "PHI-4" in chip
    assert systempage._habitability_chip({"equipment_tier": None}) == ""
    assert systempage._habitability_chip({}) == ""


def test_a_molten_or_deep_frozen_world_scores_without_overflowing():
    for kelvin in (3.0, 40.0, 1500.0, 4000.0, 20000.0):
        scored = hw.compute(_earth(surface_temperature_k=kelvin, hydrosphere="dry", p_h2o_kpa=0.0))
        assert scored["tier_temperature"] == 3
        assert 0.0 <= scored["phi_tech"] <= 1.0
