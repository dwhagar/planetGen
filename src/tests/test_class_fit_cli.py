"""The shadow scorer of the best-fit classes (GEN.1206): the report on bodies given as stored-column mappings."""

from tests import class_fit_cli as cli
from tests.test_planet_classes import CASES, body


def records():
    bodies = [(values, False) for values in CASES.values()]
    bodies.append((body(surface_temperature_k=414.0, atmospheric_pressure_pa=1000.0, ocean_fraction=0.0, hydrosphere="vapour", gravity_g=1.4), False))
    return cli.score_bodies(bodies)


def test_the_confusion_table_counts_every_body_and_agrees_on_the_core_cases():
    recs = records()
    table = cli.confusion(recs)
    assert sum(table.values()) == len(CASES) + 1
    assert all(table.get((letter, letter), 0) == 1 for letter in CASES)
    assert table[("H", "-")] == 1


def test_the_report_names_the_hole_and_the_temperate_ocean_table():
    text = cli.format_report(records(), 1, 0, 0.0)
    assert "Confusion table" in text and "Low-match list" in text
    assert "tree H, best fit -" in text
    assert "Temperate (250 to 330 K) bodies with a surface ocean" in text


def test_the_temperate_ocean_count_reads_the_hydrosphere():
    by_tree, by_fit = cli.temperate_surface_ocean(records())
    assert by_tree["M"] == 1 and by_fit["M"] == 1
