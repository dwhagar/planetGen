"""The best-fit class sweep (GEN.1206, GEN.1207): the report on bodies given as stored-column mappings."""

from tests import class_fit_cli as cli
from tests.test_planet_classes import CASES, body


def records():
    bodies = [(values, False) for values in CASES.values()]
    bodies.append((body(ocean_fraction=0.97), False))                  # a hole: a wet, mild world between the M and the O boxes
    bodies.append((body(gravity_g=0.6, mass_earth=0.3), False))        # limited by its gravity: M9 (a light world; M reaches 2 g since GEN.1306)
    return cli.score_bodies(bodies)


def test_every_case_gets_its_letter_at_100_percent_and_the_hole_and_the_limited_world_are_counted():
    recs = records()
    assert len(recs) == len(CASES) + 2
    assert all(r["match"] == 1.0 for r in recs[:len(CASES)]) and [r["fit"] for r in recs[:len(CASES)]] == list(CASES)
    assert recs[-2]["match"] == 0.0 and recs[-1]["digit"] == 9


def test_the_report_names_the_hole_the_limited_world_and_the_temperate_ocean_table():
    text = cli.format_report(records(), 1, 0, 0.0)
    assert "Share of each class code" in text and "Low-match list" in text
    assert "1 fit no class" in text and "1 are limited by another domain" in text
    assert "Temperate (250 to 330 K) bodies with a surface ocean" in text


def test_the_temperate_ocean_count_reads_the_hydrosphere():
    assert cli.temperate_surface_ocean(records())["M"] >= 1
