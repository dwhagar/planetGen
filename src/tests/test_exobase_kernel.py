"""
The compiled exobase solve (PERF.291, `planetgen.physics.exobase_kernel`): the same bits as the pure-Python solve it replaced, which is kept here as the reference.
"""

import math

import numpy as np
import pytest

from planetgen.physics import thermosphere

R_EARTH = 6.371e6


def _reference(radius_m, gravity_ms2, mean_mass, temperature_k, mesopause_k, max_radii, max_panels, bates_heights, mesopause_pa, cross_section_m2):
    """`thermosphere._exobase_solved` as it was before the kernel (CPython arithmetic)."""
    over_mass_gravity = mean_mass * gravity_ms2
    length = bates_heights * thermosphere._K * mesopause_k / over_mass_gravity
    inv_length = 1.0 / length
    limit = (max_radii - 1.0) * radius_m
    rise = temperature_k - mesopause_k
    exp, log = math.exp, math.log
    coefficient = over_mass_gravity * radius_m ** 2 / thermosphere._K                     # `1 / H = coefficient / ((R + z)^2 T)`
    log_offset = log(mesopause_pa / (thermosphere._K * mesopause_k)) + log(mesopause_k) + log(cross_section_m2)       # `f = log_offset - I - ln T - ln (1 / H)`
    panel_scale = thermosphere.EXOBASE_PANEL_SCALE_HEIGHTS
    inner = math.sqrt(0.2)                                                  # Gauss-Lobatto of four points: nodes -1, -1/sqrt(5), 1/sqrt(5), 1 with weights 1/6, 5/6, 5/6, 1/6

    z = integral = 0.0
    decay = 1.0
    temperature = mesopause_k
    inverse_height = coefficient / (radius_m * radius_m * temperature)
    value = log_offset - log(temperature) - log(inverse_height)
    height = None
    if value <= 0.0:
        height = 0.0
    else:
        bracket = None
        for _ in range(max_panels):
            rise_decay = rise * decay
            width = min(1.0 / inverse_height, length, 0.5 * (radius_m + z))
            if rise_decay > 0.0:
                width = min(width, temperature * length / rise_decay)
            z_end = min(z + panel_scale * width, limit)
            half = 0.5 * (z_end - z)
            middle = z + half
            offset = inner * half
            z_a, z_b = middle - offset, middle + offset
            t_a = temperature_k - rise * exp(-z_a * inv_length)
            t_b = temperature_k - rise * exp(-z_b * inv_length)
            decay_end = exp(-z_end * inv_length)
            t_end = temperature_k - rise * decay_end
            radius_end = radius_m + z_end
            inverse_end = coefficient / (radius_end * radius_end * t_end)
            integral_end = integral + half * ((inverse_height + inverse_end) / 6.0
                                              + (5.0 / 6.0) * (coefficient / ((radius_m + z_a) ** 2 * t_a) + coefficient / ((radius_m + z_b) ** 2 * t_b)))
            value_end = log_offset - integral_end - log(t_end) - log(inverse_end)
            if value_end <= 0.0:
                bracket = (z, integral, value, z_end, integral_end, value_end)
                break
            if z_end >= limit:
                break
            z, decay, temperature, inverse_height, integral, value = z_end, decay_end, t_end, inverse_end, integral_end, value_end
        if bracket is None:
            height = limit         # the profile does not close inside `MAX_EXOBASE_RADII` (or the panels ran out): the air is taken to reach the limit
            if z_end < limit:
                height = z_end
        else:
            low_z, low_i, low_f, high_z, high_i, high_f = bracket
            height = high_z
            for _ in range(100):
                base_z, base_i, base_f = (high_z, high_i, high_f) if abs(high_f) < abs(low_f) else (low_z, low_i, low_f)
                base_t = temperature_k - rise * exp(-base_z * inv_length)
                base_inverse = coefficient / ((radius_m + base_z) ** 2 * base_t)
                slope = 2.0 / (radius_m + base_z) - base_inverse
                x = base_z - base_f / slope if slope != 0.0 else 0.5 * (low_z + high_z)
                if not low_z < x < high_z:
                    x = 0.5 * (low_z + high_z)
                half = 0.5 * (x - base_z)
                middle = base_z + half
                offset = inner * half
                z_a, z_b = middle - offset, middle + offset
                t_a = temperature_k - rise * exp(-z_a * inv_length)
                t_b = temperature_k - rise * exp(-z_b * inv_length)
                t_x = temperature_k - rise * exp(-x * inv_length)
                inverse_x = coefficient / ((radius_m + x) ** 2 * t_x)
                integral_x = base_i + half * ((base_inverse + inverse_x) / 6.0
                                              + (5.0 / 6.0) * (coefficient / ((radius_m + z_a) ** 2 * t_a) + coefficient / ((radius_m + z_b) ** 2 * t_b)))
                value_x = log_offset - integral_x - log(t_x) - log(inverse_x)
                if value_x > 0.0:
                    low_z, low_i, low_f = x, integral_x, value_x
                else:
                    high_z, high_i, high_f = x, integral_x, value_x
                converged = abs(x - height) <= 1e-12 * x or high_z - low_z <= 1e-12 * high_z or value_x == 0.0
                height = x
                if converged:
                    break
    temperature_at = temperature_k - rise * exp(-height / length)
    scale = thermosphere._K * temperature_at / (over_mass_gravity * (radius_m / (radius_m + height)) ** 2)
    return radius_m + height, 1.0 / (cross_section_m2 * scale), scale


def _inputs(rng, count):
    """Planet-like exobase inputs: radius, gravity, mean mass, exobase and mesopause temperature."""
    for _ in range(count):
        radius = R_EARTH * float(np.exp(rng.uniform(math.log(0.2), math.log(12.0))))
        gravity = 9.8 * float(np.exp(rng.uniform(math.log(0.2), math.log(8.0))))
        mass = 1.66e-27 * float(rng.uniform(2.0, 44.0))
        mesopause = float(rng.uniform(60.0, 600.0))
        yield radius, gravity, mass, mesopause + float(np.exp(rng.uniform(math.log(10.0), math.log(6000.0)))), mesopause


def _constants():
    return (thermosphere.MAX_EXOBASE_RADII, thermosphere.MAX_EXOBASE_PANELS, thermosphere.BATES_LENGTH_SCALE_HEIGHTS, thermosphere.MESOPAUSE_PRESSURE_PA,
            thermosphere.COLLISION_CROSS_SECTION_M2)


def test_the_kernel_equals_the_python_solve_bit_for_bit_on_random_airs():
    rng = np.random.default_rng(291)
    for radius, gravity, mass, temperature, mesopause in _inputs(rng, 6000):
        want = _reference(radius, gravity, mass, temperature, mesopause, *_constants())
        assert thermosphere._exobase_solved.__wrapped__(radius, gravity, mass, temperature, mesopause, *_constants()) == want


def test_a_hot_light_air_that_never_closes_and_a_cold_one_that_closes_at_once():
    for args in ((R_EARTH, 9.8, 1.66e-27 * 2.0, 30000.0, 250.0), (R_EARTH, 9.8, 1.66e-27 * 44.0, 251.0, 250.0), (0.3 * R_EARTH, 1.0, 1.66e-27, 5000.0, 100.0)):
        assert thermosphere._exobase_solved.__wrapped__(*args, *_constants()) == _reference(*args, *_constants())


def test_the_exobase_of_earthlike_air_is_where_it_was():
    radius, density, scale = thermosphere.exobase(R_EARTH, 9.8, 28.0 * 1.66e-27, 1000.0, 200.0)
    assert 6.6e6 < radius < 7.5e6 and density > 0.0 and scale > 0.0


def test_the_inputs_must_still_be_finite_and_positive():
    with pytest.raises(ValueError):
        thermosphere._exobase_solved.__wrapped__(R_EARTH, 9.8, float("nan"), 1000.0, 200.0, *_constants())
