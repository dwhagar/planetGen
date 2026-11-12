"""UX.297, GEN.1532: the property sheet names the stored pressure of a body with no ground the reference-level pressure, and a body with a ground keeps its surface pressure."""

import pytest

from planetgen.db import store
from planetgen.describe import sheet
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw


def _pressure_label(conn, kind, row_id):
    page = sheet.sheet(conn, kind, row_id)
    labels = [row["label"] for group in page["groups"] for row in group["rows"] if row["key"] == "atmospheric_pressure_pa"]
    assert len(labels) == 1
    return labels[0]


def test_a_giant_says_reference_level_and_a_rocky_planet_says_surface(mysql_config):
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    config.MOONS = False
    config.BINARY_SYSTEM = False
    for seed in range(1, 60):
        draw.set_run_seed(seed)
        system = StarSystem(system_config=config)
        if any(b.body_type == "g" for b in system.planets) and any(b.body_type == "t" for b in system.planets):
            saved = store.save_system(system, config, config=mysql_config)
            break
    else:
        pytest.fail("no system with a giant and a rocky planet in 59 seeds")
    with store.get_connection(mysql_config) as conn:
        giant = conn.execute("SELECT id FROM planets WHERE star_system_id = ? AND body_type = 'g'", (saved,)).fetchone()["id"]
        rocky = conn.execute("SELECT id FROM planets WHERE star_system_id = ? AND body_type = 't' AND planet_type = 'terrestrial'", (saved,)).fetchone()
        assert _pressure_label(conn, "planet", giant) == "Pressure at the 1 bar reference level"
        if rocky is not None:
            assert _pressure_label(conn, "planet", rocky["id"]) == "Surface pressure"
