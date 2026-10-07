# tests/test_sector_summary.py

"""
The per-sector summary `generate.py` prints after each saved sector:
white dwarfs, giants and compact remnants under their own labels (UX.34),
and one log record per line so each keeps the debug log's prefix (OPS.9).
"""

from types import SimpleNamespace

import generate


def _entry(star_type, yerkes):
    star = SimpleNamespace(type=star_type, yerkes_class=yerkes)
    return SimpleNamespace(star_system=SimpleNamespace(star=star))


def _sector(*entries):
    return SimpleNamespace(entries=list(entries), phenomena=[], expected_system_count=lambda: 10.0)


def test_white_dwarfs_giants_and_remnants_get_their_own_labels():
    sector = _sector(
        _entry("B2D White Dwarf Star", "D"), _entry("A0D White Dwarf Star", "D"),
        _entry("A1V White Main Sequence Star", "V"), _entry("K3III Orange Giant Star", "III"),
        _entry("M2IA Red Supergiant Star", "IA"), _entry("Stellar-Mass Black Hole", "BH"),
        _entry("Neutron Star", "NS"),
    )
    lines = generate.sector_generation_summary_lines(sector, SimpleNamespace(density=1.0, num_systems=None))
    systems = lines.splitlines()[0]
    assert systems == ("    Systems: 1 A-type, 1 K-type giant, 1 M-type supergiant, 1 black hole, "
                       "1 neutron star, 2 white dwarfs")
    assert "B-type" not in systems


def test_the_summary_is_logged_one_line_per_record(monkeypatch):
    logged = []
    monkeypatch.setattr(generate.log, "normal", logged.append)
    generate._log_summary("    Systems: 1 G-type\n    Star density: actual 1.00x local, expected 1.00x local")
    assert logged == ["    Systems: 1 G-type", "    Star density: actual 1.00x local, expected 1.00x local"]
