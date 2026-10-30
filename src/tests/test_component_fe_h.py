# tests/test_component_fe_h.py

"""
GEN.662: the stars of a cluster or a satellite draw their birth [Fe/H] from their own component's distribution, store it,
and no longer draw the halo's (or read a reference 0) on the way to their lifetimes.
"""

import statistics

import pytest

from planetgen import tuning
from planetgen.galaxy import component_field, components
from planetgen.generation import bright_stars, star_population
from planetgen.physics import metallicity, stellar_evolution as evolution
from planetgen.util import draw


def _component(kind, fe_h, profile="king"):
    return components.Component(
        kind=kind, name=f"{kind}-{fe_h}", source="generated", center_pc=(0.0, 0.0, 0.0), velocity_kms=(0.0, 0.0, 0.0),
        epoch_unix=0.0, mass_solar=1.0e5, half_mass_radius_pc=5.0, profile=profile, age_gy=11.0, metallicity_fe_h=fe_h,
        seed=1, tidal_radius_pc=50.0, king_w0=6 if profile == "king" else None)


def test_a_component_supplies_its_own_value_with_the_scatter_of_its_family():
    cluster = _component("cluster", -1.6)
    assert component_field.fe_h_peaks_of(cluster) == ((1.0, -1.6, tuning.COMPONENT_FE_H_SIGMA["cluster"]),)
    satellite = _component("satellite", -0.7, profile="exponential")
    assert component_field.fe_h_peaks_of(satellite) == ((1.0, -0.7, tuning.COMPONENT_FE_H_SIGMA["satellite"]),)


def test_the_stars_of_a_component_centre_on_its_value_and_are_stored():
    peaks = component_field.fe_h_peaks_of(_component("cluster", -1.6))
    stars = star_population.sample_bright_stars(
        300, 50.0, "cluster", draw.Stream("gen662"),
        fe_h_source=bright_stars.spot_fe_h("cluster", [((0, 0, 0),)] * 300, draw.Stream("gen662-fe"), peaks))
    values = [star["metallicity_fe_h"] for star in stars]
    assert statistics.mean(values) == pytest.approx(-1.6, abs=0.1)
    assert all(bright_stars.star_fe_h(star) == star["metallicity_fe_h"] for star in stars)


def test_a_cluster_or_satellite_star_with_no_component_takes_the_generic_distribution_not_the_halos():
    stream = draw.Stream("generic")
    cluster = [evolution.birth_fe_h(stream, "cluster") for _ in range(400)]
    satellite = [evolution.birth_fe_h(stream, "satellite_young") for _ in range(400)]
    assert statistics.mean(cluster) == pytest.approx(
        sum(w * m for w, m, _s in tuning.METALLICITY_CLUSTER_PEAKS) / sum(w for w, _m, _s in tuning.METALLICITY_CLUSTER_PEAKS),
        abs=0.2)
    assert statistics.mean(satellite) == pytest.approx(tuning.METALLICITY_SATELLITE_PEAKS[0][1], abs=0.15)
    assert metallicity.generic_peaks("halo") == ()


def test_the_stars_of_two_components_follow_their_own_values():
    poor = [evolution.birth_fe_h(draw.Stream(f"p{k}"), "cluster", component_field.fe_h_peaks_of(_component("cluster", -2.0)))
            for k in range(200)]
    rich = [evolution.birth_fe_h(draw.Stream(f"r{k}"), "cluster", component_field.fe_h_peaks_of(_component("cluster", -0.5)))
            for k in range(200)]
    assert statistics.mean(rich) - statistics.mean(poor) == pytest.approx(1.5, abs=0.1)


def test_a_components_stars_are_kept_against_the_best_composition_they_can_have():
    """A cluster far from the best bin of the tables still gets its stars: the source is stamped with its support."""
    peaks = component_field.fe_h_peaks_of(_component("cluster", -2.1))
    source = bright_stars.spot_fe_h("cluster", [((0, 0, 0),)], draw.Stream("support"), peaks)
    assert source.bins and star_population.fe_h_bin(-2.1) in source.bins
    assert star_population.fe_h_bin(0.0) not in source.bins
    stars = star_population.sample_bright_stars(5, 500.0, "cluster", draw.Stream("far"), fe_h_source=source)
    assert len(stars) == 5 and all(star["metallicity_fe_h"] < -1.0 for star in stars)
