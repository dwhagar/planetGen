# tests/test_taxon_kernel.py

"""GEN.473: the taxa of a belt from the disk's ice lines and the early heating of its bodies (`physics/taxon_kernel.py`)."""

import math

import pytest

from planetgen.physics import condensate_shares, small_body_composition as sb, taxon_kernel as tk

SUN_W = 3.828e26
SOLAR_SYSTEM_RATIO = 5.23e-5      # a test point ([R2128] via [R2768]), never an input

FAMILY_OF = {"E": "dry", "S": "dry", "A": "dry", "L": "dry", "K": "dry", "V": "dry", "M": "dry", "C": "hydrated", "B": "hydrated",
             "P": "primitive", "D": "ice-rich"}


def test_the_families_change_at_the_ice_lines_of_lodders():
    assert tk.WATER_LINE_K == 182.0 and tk.AMMONIA_LINE_K == 131.0
    assert tk.family_shares(260.0)["dry"] > 0.9
    at_water = tk.family_shares(tk.WATER_LINE_K, implant=0.0)
    assert at_water["dry"] == pytest.approx(0.5)
    at_ammonia = tk.family_shares(tk.AMMONIA_LINE_K, implant=0.0)
    assert at_ammonia["ice-rich"] == pytest.approx(0.5, abs=0.01)
    assert tk.family_shares(100.0)["ice-rich"] > 0.9


def test_the_shares_always_add_to_one_and_the_kernel_too():
    for temperature in (30.0, 100.0, 131.0, 150.0, 166.0, 182.0, 250.0, 900.0):
        shares = tk.family_shares(temperature)
        assert sum(shares.values()) == pytest.approx(1.0)
        assert all(value >= 0.0 for value in shares.values())
        assert sum(tk.kernel(temperature, 2e-5).values()) == pytest.approx(1.0)


def test_the_boundaries_move_with_the_stars_luminosity():
    """The temperature at a distance is the star's own, so a bright star's belt is the same one pushed out by the square root of the luminosity."""
    sun_line = condensate_shares.ice_line_au(SUN_W)
    bright = condensate_shares.ice_line_au(4.0 * SUN_W)
    assert bright == pytest.approx(2.0 * sun_line)
    inner_sun = tk.kernel_at_distance(0.8 * sun_line, SUN_W, 1e-6)
    inner_bright = tk.kernel_at_distance(0.8 * bright, 4.0 * SUN_W, 1e-6)
    for taxon in tk.TAXA:
        assert inner_sun[taxon] == pytest.approx(inner_bright[taxon])


def test_the_melting_ratio_is_a_few_times_the_galactic_mean_and_a_third_of_the_solar_systems():
    assert 5e-6 < tk.melting_ratio() < 3e-5
    assert tk.melted_share(0.5 * tk.melting_ratio()) == 0.0
    assert tk.melted_share(tk.melting_ratio()) == 0.0
    assert 0.3 < tk.melted_share(SOLAR_SYSTEM_RATIO) < 0.8
    assert tk.melted_share(1e-3) == 1.0
    values = [tk.melted_share(r) for r in (1e-5, 2e-5, 4e-5, 8e-5, 1.6e-4)]
    assert values == sorted(values)


def test_a_belt_round_a_star_born_without_aluminium_26_has_no_igneous_taxa():
    cold = tk.kernel_at_distance(2.0, SUN_W, 1e-7)
    assert cold["V"] == cold["A"] == cold["M"] == 0.0
    assert cold["S"] > 0.5
    hot = tk.kernel_at_distance(2.0, SUN_W, 2e-4)
    assert hot["S"] == hot["K"] == hot["L"] == 0.0
    assert hot["A"] > hot["M"] > hot["V"] > 0.0


def test_the_solar_systems_zones_come_back_as_one_realisation_of_the_derived_kernel():
    """Table 6 of [R1065], summed by family, against the derived families at each zone's temperature (rms of the shares under 0.2)."""
    errors = []
    for name, low, high, shares in sb.ZONES:
        x = math.sqrt(low * high)
        temperature = tk.WATER_LINE_K / math.sqrt(x)
        total = float(sum(shares.values()))
        observed = {family: sum(v for taxon, v in shares.items() if FAMILY_OF[taxon] == family) / total for family in tk.FAMILIES}
        modelled = tk.family_shares(temperature)
        errors += [observed[family] - modelled[family] for family in tk.FAMILIES]
    assert math.sqrt(sum(e * e for e in errors) / len(errors)) < 0.2


def test_the_trojans_sit_at_the_ammonia_line_and_the_belt_changes_at_the_snow_line():
    """A prediction and not a fit: the ammonia line is at 1.93 snow-line distances and the Trojans (5.2 AU, 1.9) are ice-rich D."""
    x_ammonia = (tk.WATER_LINE_K / tk.AMMONIA_LINE_K) ** 2
    assert x_ammonia == pytest.approx(1.93, abs=0.01)
    assert tk.family_shares(tk.WATER_LINE_K / math.sqrt(1.93))["ice-rich"] > 0.35
    assert tk.family_shares(tk.WATER_LINE_K / math.sqrt(0.8))["dry"] > 0.7
    assert tk.family_shares(tk.WATER_LINE_K / math.sqrt(1.1))["hydrated"] > 0.3


def test_bad_inputs_raise():
    with pytest.raises(ValueError):
        tk.family_shares(0.0)
