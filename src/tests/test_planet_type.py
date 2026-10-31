"""The body type (GEN.1219): terrestrial, transitional or gas by the envelope share of the mass."""

import pytest

from planetgen.physics import planet_type as pt
from tests import planet_type_cli as cli


@pytest.mark.parametrize("envelope, kind", [(0.0, "terrestrial"), (0.0049, "terrestrial"), (0.005, "transitional"), (0.15, "transitional"),
                                            (0.4999, "transitional"), (0.5, "gas"), (0.95, "gas")])
def test_the_cuts_are_clean(envelope, kind):
    assert pt.type_of(envelope) == (kind, True)


def test_a_body_with_no_stored_envelope_is_typed_by_body_type_and_not_known():
    assert pt.type_of(None, "g") == ("gas", False) and pt.type_of(None, "t") == ("terrestrial", False)


def test_planet_type_reads_an_object_or_a_mapping():
    class Body:
        envelope_mass_fraction, body_type = 0.12, "g"
    assert pt.planet_type(Body()) == "transitional"
    assert pt.planet_type({"envelope_mass_fraction": 0.9, "body_type": "g"}) == "gas"


def test_the_gas_cut_is_the_giant_line_of_the_bulk_structure():
    from planetgen.physics import bulk_structure
    assert pt.GAS_MIN_ENVELOPE == bulk_structure.GIANT_ENVELOPE_FRACTION


def test_the_report_counts_every_body():
    rows = [{"type": "terrestrial", "known": True, "envelope": 0.0, "body_type": "t", "class": "C", "moon": False, "radius": 1.0, "water": 0.0, "mass": 1.0},
            {"type": "gas", "known": True, "envelope": 0.9, "body_type": "g", "class": "J", "moon": False, "radius": 10.0, "water": 0.0, "mass": 300.0}]
    text = cli.format_report(rows, 1, 0, 0.0)
    assert "| terrestrial | 1 | 0 | 0 | 1 |" in text and "| gas | 0 | 1 | 0 | 1 |" in text
