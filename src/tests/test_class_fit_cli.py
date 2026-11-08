"""The best-fit class sweep (GEN.1206, GEN.1207): the report on bodies given as stored-column mappings."""

from tests import class_fit_cli as cli
from tests.test_planet_classes import CASES, body


def records():
    bodies = [(values, False) for values in CASES.values()]
    bodies.append((body(surface_temperature_k=310.0, ocean_fraction=0.5), False))   # a hole: a half-ocean world past the M box (GEN.1342 closed the wet, mild gap to O)
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


def test_a_table_compared_with_itself_changes_no_class_and_one_with_a_moved_edge_changes_some():
    """GEN.1321: `--compare` counts the bodies whose class code differs under another fit table."""
    import json
    from planetgen.physics import class_fit
    recs = records()
    same = open(class_fit.DATA_FILE, encoding="utf-8").read()
    changed, total, pairs = cli.changes_against(recs, same)
    assert changed == 0 and total >= len(CASES) and not pairs
    document = json.loads(same)
    for cls in document["classes"]:
        for condition in cls["when"]:
            if cls["letter"] == "M" and condition["quantity"] == "o2_partial_pa" and "min" in condition:
                condition["min"] = 1.0e9                      # M no longer fits any body: they take other classes
    changed, total, pairs = cli.changes_against(recs, json.dumps(document))
    assert changed >= 1 and sum(pairs.values()) == changed
    assert class_fit.CLASSES is not None and "o2_partial_pa" in {c["quantity"] for k in class_fit.CLASSES if k.letter == "M" for c in k.when}      # the table is put back
    assert "## Class codes against the fit table at X" in cli.format_changes("X", changed, total, pairs)
