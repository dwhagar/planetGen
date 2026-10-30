"""
The Pecaut and Mamajek dwarf spectral-type sequence (GEN.943, `datasets/mamajek_sequence.py`) and what is built on it.
"""
import pytest

from planetgen import tuning
from planetgen.datasets import fetch, mamajek_sequence
from planetgen.physics import constants, stellar_evolution

TABLE = """#SpT   Teff  logT   BCv    logL   Mbol R_Rsun   Mv    B-V
O3V    44900 4.652 -4.01   5.82  -9.81 13.43   -5.80 -0.330
G2V    5770 3.761 -0.085  0.00 4.82 1.0 4.8 0.65
G2.5V  5750 3.759 -0.09  0.00 4.82 1.0 4.8 0.65
M9V    2380 3.377 -5.86  -3.9 14.6 0.1 15 2.0
#SpT   Teff  logT  BCv
B0V    31400 4.497 -2.99 4.0 0 0 0 0
"""


def test_the_parser_reads_whole_number_types_of_the_first_table_only():
    with pytest.raises(ValueError, match="not the 67"):
        mamajek_sequence.parse(TABLE)


def test_the_real_file_holds_sixty_seven_types_from_o3_to_m9():
    rows = [row for row in mamajek_sequence._load() if row.log_luminosity is not None]
    assert len(rows) == 67
    assert (rows[0].letter, rows[0].subclass, rows[0].teff_k) == ("O", 3, 44900.0)
    assert (rows[-1].letter, rows[-1].subclass) == ("M", 9)
    assert all(a.teff_k > b.teff_k for a, b in zip(rows, rows[1:]))


def test_the_hot_o_types_are_the_stand_ins():
    hot = [row for row in mamajek_sequence._load() if row.log_luminosity is None]
    assert {row.subclass: row.teff_k for row in hot} == tuning.HOT_O_SUBCLASS_TEMPERATURES_K


def test_every_type_reads_back_as_itself():
    for row in mamajek_sequence._load():
        assert mamajek_sequence.spectral_letter_and_subclass(row.teff_k) == (row.letter, row.subclass)
        assert mamajek_sequence.temperature_k(row.letter, row.subclass) == row.teff_k


def test_the_letter_ranges_are_the_map_bands_of_the_constants():
    for letter, (low, high) in constants.TEMP_RANGES.items():
        found = mamajek_sequence.letter_range_k(letter)
        assert found == pytest.approx((low, high), abs=1.0)


def test_the_sun_is_not_a_subdwarf_and_a_metal_poor_star_is():
    sun = {"yerkes_class": "V", "luminosity_sol": 1.0, "temperature_k": 5772.0}
    assert stellar_evolution.yerkes_class_of(sun) == "V"
    poor = dict(sun, luminosity_sol=0.3)
    assert stellar_evolution.yerkes_class_of(poor) == "VI"
    assert stellar_evolution.yerkes_class_of(dict(poor, yerkes_class="III")) == "III"


@pytest.mark.real_datasets
def test_a_missing_file_is_refused(tmp_path, monkeypatch):
    from tests import datasets_support as support
    support.reset_caches(tmp_path / "cache", monkeypatch)
    with pytest.raises(fetch.DataMissing):
        mamajek_sequence._load()
