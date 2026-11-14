"""
The stand-in ranking by share of the fudge (GEN.1316, scripts/build_fudge_ranking.py) and the softness of the class edges on the class pages (class_fit.softness_words).
"""

import json
import pathlib
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import build_fudge_ranking as ranking  # noqa: E402
import build_fudge_table as table_builder  # noqa: E402
import stand_in_sensitivity  # noqa: E402

from planetgen.physics import class_fit  # noqa: E402

MATRIX = json.loads(pathlib.Path(table_builder.MATRIX).read_text(encoding="utf-8"))
TABLE = json.loads(pathlib.Path(table_builder.OUT).read_text(encoding="utf-8"))
ROWS = stand_in_sensitivity.read_rows()


def test_the_page_is_the_table_as_built():
    assert pathlib.Path(ranking.OUT).read_text(encoding="utf-8") == ranking.render(MATRIX, TABLE, ROWS), "run python scripts/build_fudge_ranking.py"


def test_the_shares_of_a_quantity_add_to_one():
    for quantity, (sigma, shared) in ranking.shares(MATRIX, TABLE, ROWS).items():
        assert sigma == pytest.approx(TABLE["quantities"][quantity]["sigma"], abs=1.0e-6), quantity
        assert sum(shared.values()) == pytest.approx(1.0, abs=1.0e-9), quantity


def test_the_shares_agree_with_the_five_the_table_names():
    found = ranking.shares(MATRIX, TABLE, ROWS)
    for quantity, entry in TABLE["quantities"].items():
        for contributor in entry["contributors"]:
            assert found[quantity][1][contributor["fill_id"]] == pytest.approx(contributor["share"], abs=1.0e-4), (quantity, contributor["fill_id"])


def test_a_quantity_with_no_measured_row_has_no_shares():
    found = ranking.shares(MATRIX, TABLE, ROWS)
    for quantity, entry in TABLE["quantities"].items():
        assert (quantity in found) is (entry["basis"] == "measured")


@pytest.mark.parametrize("tier,tweak,score", [("T1", "L0", 0.0), ("T9", "L4", 12.0), ("T6", "L1", 5.5), ("T4", "L3", 6.0)])
def test_the_unreliability_score_is_the_ranking_reports(tier, tweak, score):
    assert ranking.unreliability(tier, tweak) == score


def test_the_page_says_what_the_fudge_is():
    text = pathlib.Path(ranking.OUT).read_text(encoding="utf-8")
    assert "the model's uncertainty, not the world's variety" in text


def test_every_planet_class_edge_has_a_softness_phrase():
    for cls in class_fit.CLASSES:
        edged = [c for c in cls.when if "equals" not in c]
        assert len(class_fit.softness_words(cls)) == len(edged), cls.letter


def test_a_softness_phrase_gives_the_width_in_the_unit_of_the_quantity():
    letter_m = next(c for c in class_fit.CLASSES if c.letter == "M")
    words = class_fit.softness_words(letter_m)
    assert "mean surface temperature at least 275 K: soft by 0.925 K" in words
    assert any("surface gravity at least 0.75 g: soft by a factor of" in line for line in words)
    assert any("percentage points" in line for line in words)


def test_the_class_pages_carry_the_softness_and_the_note():
    from planetgen.web.lib import classref

    entry = classref.catalog()["planet"]
    assert any("not the variety of the worlds" in note for note in entry["notes"])
    for code, item in entry["classes"].items():
        assert any(label.startswith("Soft edges") for label, _ in item["facts"]), code
