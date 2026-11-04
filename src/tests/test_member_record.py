# tests/test_member_record.py

"""
GEN.1009: the member record (`planetgen/generation/member_record.py`), one frozen evolution of a member that the tree and the
build both read (docs/design/systems/overlap-free-systems.md section 4.1).
"""

import dataclasses

import pytest

from planetgen.generation import member_record
from planetgen.physics import alpha_tracks, constants, post_main_sequence
from planetgen.physics import unit_registry as units

_km_to_au = units.converter("km", "au")

GY = 1.0e9


def test_the_record_is_a_function_of_seed_mass_age_and_fe_h():
    first = member_record.make(11, 1.4, 3.0 * GY, 0.1)
    assert first == member_record.make(11, 1.4, 3.0 * GY, 0.1)
    assert hash(first) == hash(member_record.make(11, 1.4, 3.0 * GY, 0.1))


def test_the_record_is_frozen():
    record = member_record.make(11, 1.0, 4.6 * GY, 0.0)
    with pytest.raises(dataclasses.FrozenInstanceError):
        record.radius_au = 1.0


def test_a_horizontal_branch_star_keeps_the_colour_its_seed_drew():
    """The colour and luminosity of a horizontal-branch star are drawn; the stream is the member's own, so a second ask agrees and
    different seeds are different draws of the same star."""
    track = alpha_tracks.track_fe_h(-1.5, 12.0)    # GEN.1024: the track a star this old reads (the horizontal branch is past the gate)
    window = post_main_sequence.timeline(0.85, track).window(post_main_sequence.HORIZONTAL_BRANCH)
    age_yr = 0.5 * (window.start_gy + window.end_gy) * GY
    drawn = [member_record.make(seed, 0.85, age_yr, -1.5) for seed in range(30)]
    assert drawn == [member_record.make(seed, 0.85, age_yr, -1.5) for seed in range(30)]
    assert {record.state for record in drawn} == {"III"}
    assert len({record.teff_k for record in drawn}) > 1


@pytest.mark.parametrize("mass, age_gy, kind", [(1.0, 4.6, "main-sequence"), (1.0, 11.0, "giant"), (2.0, 1.5, "white-dwarf"),
                                               (1.0, 0.002, "pre-main-sequence"), (0.05, 1.0, "brown-dwarf")])
def test_the_kind_follows_the_evolution(mass, age_gy, kind):
    assert member_record.make(3, mass, age_gy * GY, 0.0).kind == kind


def test_the_radius_is_the_one_the_star_is_made_with():
    """`radius_au` is `star_params`' radius (from the temperature rounded to a hundred kelvin), so a star made from the record has
    exactly that radius."""
    for mass, age_gy in ((1.0, 4.6), (1.5, 3.0), (2.0, 1.5), (0.3, 5.0), (3.0, 0.4)):
        record = member_record.make(5, mass, age_gy * GY, 0.0)
        assert record.radius_au == pytest.approx(_km_to_au(record.star_params()["radius_km"]), rel=1.0e-12)


def test_a_white_dwarf_remembers_its_giant_tip():
    record = member_record.make(9, 2.0, 1.5 * GY, 0.0)
    assert record.kind == "white-dwarf"
    assert record.radius_max_au > 1000.0 * record.radius_au
    assert record.mass_max_solar == pytest.approx(2.0)
    assert record.mass_solar < record.mass_max_solar


def test_a_giant_on_the_branch_after_the_first_tip_keeps_the_tip():
    record = member_record.make(9, 1.5, 3.0 * GY, 0.1)
    assert record.kind == "giant"
    assert record.radius_max_au > record.radius_au


def test_a_star_on_the_main_sequence_has_no_larger_past():
    for mass, age_gy in ((1.0, 4.6), (0.5, 2.0), (1.0, 0.002)):
        record = member_record.make(9, mass, age_gy * GY, 0.0)
        assert record.radius_max_au == record.radius_au
        assert record.mass_max_solar == record.initial_mass_solar


def test_the_radius_max_helper_is_the_tip_of_the_finished_windows():
    line = post_main_sequence.timeline(1.0, 0.0)
    assert line.radius_max_sol(line.main_sequence_gy * 0.5) == 0.0
    after = line.radius_max_sol(line.end_gy * 2.0)
    assert after > 100.0                                  # a one solar-mass star's red-giant tip is a couple of hundred radii
    assert line.radius_max_sol(line.end_gy * 3.0) == after
    midway = line.radius_max_sol(line.windows[1].start_gy + 0.5 * (line.windows[1].end_gy - line.windows[1].start_gy))
    assert 0.0 < midway < after


def test_a_collapsed_birth_mass_is_lightened_or_raises():
    lightened = member_record.make(2, 30.0, 12.0 * GY, 0.0)
    assert lightened.initial_mass_solar < 30.0
    assert lightened.mass_max_solar == lightened.initial_mass_solar
    with pytest.raises(ValueError):
        member_record.make(2, 30.0, 12.0 * GY, 0.0, lighten=False)


def test_a_brown_dwarf_has_no_star_params():
    record = member_record.make(4, 0.05, 1.0 * GY, 0.0)
    assert record.star_params() is None
    assert record.mass_solar == record.initial_mass_solar == record.mass_max_solar


def test_a_compact_remnant_keeps_its_progenitors_mass():
    record = member_record.compact("neutron-star", 8, 1.4, 12.0, 1.0e-7, 1.0e8, 0.0)
    assert (record.kind, record.state, record.params) == ("neutron-star", None, None)
    assert record.mass_max_solar == 12.0
    assert record.radius_max_au == record.radius_au == 1.0e-7
