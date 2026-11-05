"""
Rogue planet classes (GEN.8, GEN.687): a rogue planet's class is derived from its mass, radius and the surface its own heat gives
it (`rogue.rogue_class`, the same `planet_classes.classify` as any body), and the class is stored (schema v47) and shown on the
phenomenon page, the sector's Contents and the Sector Map.
"""


from planetgen.web.app import create_app
from planetgen.api.config import Config

from planetgen import tuning
from planetgen.generation import phenomena_plausibility as pp
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import RoguePlanet

from planetgen.web.maps import starmap  # noqa: E402
from planetgen.web import sector_page, system_pages  # noqa: E402
import pytest
from tests import render_support


class _Config(Config):
    TESTING = True
    SECRET_KEY = "test"


def test_no_class_depends_on_a_rogue_flag_and_none_has_life():
    for code, data in tuning.PLANET_CLASSES.items():
        assert "r" not in data, code
        assert not data.get("life_class"), code


def test_brown_dwarfs_have_no_planet_class():
    assert RoguePlanet(SystemConfig()).planet_class is None


def test_plausibility_flags_a_class_of_the_other_type():
    record = {"phenomenon_type": "rogue-planet", "mass_kg": 5e28, "radius_km": 70000.0, "planet_type": "g",
              "planet_class": "M", "galactic_orbital_speed_kms": 205.7, "galactic_orbital_period_gy": 0.236,
              "galactic_orbital_phase_deg": 180.0, "galactic_min_update_interval_years": 1e-9}
    assert any("planet_class" in issue for issue in pp.check_hard_invariants(record))
    record["planet_class"] = None
    assert not any("planet_class" in issue for issue in pp.check_hard_invariants(record))


def test_phenomenon_page_contents_and_sector_map_show_the_class():
    with create_app(_Config).test_request_context("/"):
        fields = {label: (text, url) for label, text, url in
                  system_pages.phenomenon_fields("rogue_planet", {"planet_class": "C", "planet_type": "t"})}
        assert fields["Planet Class"][0] == "C"
        assert fields["Planet Class"][1]

        rogue = {"type": "rogue_planet", "id": 1, "name": "Wanderer", "descriptor": "terrestrial", "class": "C",
                 "radius_ly": 0, "distance_ly": 1.0}
        group = sector_page._rogue_group_row([rogue, {**rogue, "id": 2, "class": None, "descriptor": "brown dwarf"}], lambda system_id: f"/s/{system_id}")
        assert [m["details"] for m in group["members"]] == ["Class C, Terrestrial", "Brown dwarf"]

    data = starmap._cloud_data(lambda *a, **k: "/x", rogue, 0, 0, 0, 1)
    assert data["typeLabel"] == "Rogue Planet (Class C, terrestrial)"
    data = starmap._cloud_data(lambda *a, **k: "/x", {**rogue, "class": None}, 0, 0, 0, 1)
    assert data["typeLabel"] == "Rogue Planet (Terrestrial)"
