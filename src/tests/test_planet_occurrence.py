# tests/test_planet_occurrence.py

"""
GEN.728 and TEST.395: the observed planet occurrence figures (`planetgen.physics.planet_occurrence`): the file and its refusals, the
host groups, the comparison of a simulated figure with a band, and the multiplicity chi-square.
"""

import json

import pytest

from planetgen.physics import planet_occurrence as po
from planetgen.physics.elements import DataError


def written(tmp_path, mutate=None):
    with open(po.DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    if mutate:
        mutate(document)
    path = tmp_path / "occurrence.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    return str(path)


def refused(tmp_path, match, mutate):
    with pytest.raises(DataError, match=match):
        po.load(written(tmp_path, mutate))


def test_every_band_names_its_paper_its_label_and_a_known_host():
    assert len(po.BANDS) >= 30
    for band in po.BANDS.values():
        assert band.host in (*po.HOSTS, "all") and band.label in po.LABELS and "[R" in band.source
        assert band.tier in ("T1", "T2") and band.unit and band.definition
        assert band.minus >= 0.0 and band.plus >= 0.0 and band.value > 0.0
    assert {b.label for b in po.BANDS.values()} == {"READ", "2ND-HAND"}
    assert all(b.tier == "T2" for b in po.BANDS.values() if b.label == "2ND-HAND")


def test_a_figure_with_no_printed_uncertainty_gets_the_default_relative_width():
    band = po.BANDS["inner_fraction.M"]
    assert not band.printed and band.minus == band.plus == pytest.approx(0.75 * 0.25)
    sun_like = po.BANDS["inner_fraction.G"]
    assert sun_like.printed and (sun_like.minus, sun_like.plus) == (0.03, 0.03)
    assert sun_like.interval(2.0) == pytest.approx((0.24, 0.36))


def test_the_bands_hold_the_figures_the_papers_print():
    assert po.BANDS["inner_fraction.G"].value == 0.30 and po.BANDS["inner_multiplicity.G"].value == 3.0
    assert po.BANDS["hot_jupiter_rate.G"].value == pytest.approx(0.0062)
    assert po.BANDS["hot_jupiter_rate.M"].value == pytest.approx(0.0027)
    assert po.BANDS["giants_2_8au"].value == pytest.approx(0.141) and po.BANDS["giants_2_8au"].plus == pytest.approx(0.020)
    assert po.BANDS["spacing_k_ge_10"].value == 0.93 and po.BANDS["spacing_k_mode"].value == 20.0
    assert po.BANDS["valley_radius_10d"].label == "2ND-HAND" and po.BANDS["valley_radius_10d"].value == 1.9


def test_the_transit_multiplicity_is_the_published_count_and_adds_up():
    multiplicity = po.TRANSIT_MULTIPLICITY
    assert multiplicity.counts == {1: 432, 2: 99, 3: 42, 4: 11, 5: 3, 6: 2}
    assert sum(multiplicity.counts.values()) == multiplicity.systems == 589
    assert multiplicity.stars_surveyed == 30759
    assert sum(multiplicity.shares().values()) == pytest.approx(1.0)


def test_a_host_belongs_to_the_group_of_its_mass():
    assert [po.host_group(m) for m in (0.2, 0.7, 1.0, 1.3, 3.0)] == ["M", "K", "G", "F", None]
    assert [b.host for b in po.bands_of("hot_jupiter_rate", "M")] == ["M"]
    assert {b.id for b in po.bands_of("hot_jupiter_rate", "G")} == {"hot_jupiter_rate.G", "hot_jupiter_rate.G.rv"}
    assert {b.id for b in po.bands_of("pair_spacing_mode_hill_radii", "G")} == {"spacing_k_mode"}


def test_a_figure_is_held_to_the_published_interval_and_its_own_sampling_error():
    band = po.BANDS["inner_fraction.G"]
    assert po.compare(0.30, band).z == 0.0
    near = po.compare(0.35, band)
    assert near.z == pytest.approx(0.05 / 0.03) and near.ok
    far = po.compare(0.40, band)
    assert far.z == pytest.approx(0.10 / 0.03) and not far.ok
    # 100 systems cannot tell 0.40 from 0.30 at two sigma once their own sampling error is counted
    assert po.compare(0.40, band, sampling_sigma=po.share_sigma(0.40, 100)).ok
    below = po.compare(0.20, band)
    assert below.z < 0 and not below.ok


def test_the_sampling_error_of_a_share_never_vanishes():
    assert po.share_sigma(0.5, 100) == pytest.approx(0.05)
    assert po.share_sigma(0.0, 100) > 0.0 and po.share_sigma(1.0, 400) > 0.0


def test_the_observed_multiplicity_matches_itself_and_a_different_one_does_not():
    chi2, dof = po.multiplicity_chi2(dict(po.TRANSIT_MULTIPLICITY.counts))
    assert chi2 < 1e-9 and dof == 4   # bin 6 pooled with 5: it expects fewer than 5
    flat = {1: 150, 2: 150, 3: 150, 4: 100, 5: 30, 6: 9}
    chi2, dof = po.multiplicity_chi2(flat)
    assert chi2 > 50
    scaled = {k: 4 * v for k, v in po.TRANSIT_MULTIPLICITY.counts.items()}
    assert po.multiplicity_chi2(scaled)[0] == pytest.approx(0.0)
    with pytest.raises(ValueError):
        po.multiplicity_chi2({})


def test_a_system_with_more_transiting_planets_than_the_survey_saw_joins_the_last_bin():
    counts = dict(po.TRANSIT_MULTIPLICITY.counts)
    counts[8] = 0
    assert po.multiplicity_chi2(counts)[0] == pytest.approx(0.0)
    counts[8] = 3
    assert po.multiplicity_chi2(counts)[0] > 0.0


# --- the refusals ----------------------------------------------------------------------------------------------------------


def test_a_missing_or_unreadable_file_is_refused(tmp_path):
    with pytest.raises(DataError, match="cannot read"):
        po.load(str(tmp_path / "none.json"))


def test_a_repeated_band_id_is_refused(tmp_path):
    refused(tmp_path, "repeated", lambda d: d["bands"].append(dict(d["bands"][0])))


def test_a_band_with_an_unknown_host_or_label_is_refused(tmp_path):
    def bad_host(d):
        d["bands"][0]["host"] = "Z"

    def bad_label(d):
        d["bands"][0]["label"] = "RECALLED"

    refused(tmp_path, "not valid", bad_host)
    refused(tmp_path, "not valid", bad_label)


def test_a_second_hand_figure_must_name_the_paper_that_quotes_it(tmp_path):
    def drop(d):
        for band in d["bands"]:
            if band["label"] == "2ND-HAND":
                band.pop("citing")

    refused(tmp_path, "citing", drop)


def test_a_band_without_a_source_or_with_a_nonfinite_value_is_refused(tmp_path):
    def no_source(d):
        d["bands"][0]["source"] = ""

    def no_value(d):
        d["bands"][0]["value"] = None

    refused(tmp_path, "source is missing", no_source)
    refused(tmp_path, "not a finite number", no_value)


def test_a_multiplicity_that_does_not_add_up_is_refused(tmp_path):
    def wrong(d):
        d["transit_multiplicity"]["systems_with_a_transiting_planet"] = 600

    refused(tmp_path, "do not add", wrong)


def test_a_host_range_that_does_not_increase_is_refused(tmp_path):
    def wrong(d):
        d["hosts"]["G"]["high"] = 0.5

    refused(tmp_path, "not increasing", wrong)
