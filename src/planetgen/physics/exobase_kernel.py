# planetgen/physics/exobase_kernel.py

"""
The compiled exobase solve (PERF.291)
=====================================

`thermosphere.exobase` for one air: the arithmetic of `_exobase_solved`, line for line, compiled. It is called 17,000 times in the 700 loss histories of one wide
binary's planets, and the pure-Python solve was a fifth of that time.

**Same bits.** The compiled function must give what CPython gives, to the last bit (`tests/test_exobase_kernel.py` holds `==` on 20,000 real calls and on random ones):
`exp`, `log` and `sqrt` already agree; a square written `x ** 2` does not, because CPython's `pow(x, 2.0)` is the library's and the compiler folds a constant exponent
to `x * x`, which differs in the last bit for about one value in a thousand. The squares here go through `math.pow(x, two)` with `two` passed in (2.0), which the compiler
cannot fold. Like `ktable_kernel`, the module is imported when the solve is called, never by the web app (numba's cache needs a writable code tree).
"""

import math

from numba import njit


@njit(cache=True)
def solve(radius_m, gravity_ms2, mean_mass, temperature_k, mesopause_k, max_radii, max_panels, bates_heights, mesopause_pa, cross_section_m2,
          panel_scale_heights, boltzmann, two):
    """`(radius of the exobase, number density there, scale height there)` as `thermosphere._exobase_solved` (the same arguments, then the panel scale in scale heights,
    the Boltzmann constant, and 2.0). The inputs are finite and positive (checked by the caller)."""
    _K = boltzmann

    over_mass_gravity = mean_mass * gravity_ms2
    length = bates_heights * _K * mesopause_k / over_mass_gravity
    inv_length = 1.0 / length
    limit = (max_radii - 1.0) * radius_m
    rise = temperature_k - mesopause_k
    coefficient = over_mass_gravity * math.pow(radius_m, two) / _K                     # `1 / H = coefficient / ((R + z)^2 T)`
    log_offset = math.log(mesopause_pa / (_K * mesopause_k)) + math.log(mesopause_k) + math.log(cross_section_m2)       # `f = log_offset - I - ln T - ln (1 / H)`
    panel_scale = panel_scale_heights
    inner = math.sqrt(0.2)                                                  # Gauss-Lobatto of four points: nodes -1, -1/sqrt(5), 1/sqrt(5), 1 with weights 1/6, 5/6, 5/6, 1/6

    z = integral = 0.0
    decay = 1.0
    temperature = mesopause_k
    inverse_height = coefficient / (radius_m * radius_m * temperature)
    value = log_offset - math.log(temperature) - math.log(inverse_height)
    height = 0.0
    if value <= 0.0:
        height = 0.0
    else:
        has_bracket = False
        low_z = low_i = low_f = high_z = high_i = high_f = 0.0
        z_end = 0.0
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
            t_a = temperature_k - rise * math.exp(-z_a * inv_length)
            t_b = temperature_k - rise * math.exp(-z_b * inv_length)
            decay_end = math.exp(-z_end * inv_length)
            t_end = temperature_k - rise * decay_end
            radius_end = radius_m + z_end
            inverse_end = coefficient / (radius_end * radius_end * t_end)
            integral_end = integral + half * ((inverse_height + inverse_end) / 6.0
                                              + (5.0 / 6.0) * (coefficient / (math.pow(radius_m + z_a, two) * t_a) + coefficient / (math.pow(radius_m + z_b, two) * t_b)))
            value_end = log_offset - integral_end - math.log(t_end) - math.log(inverse_end)
            if value_end <= 0.0:
                has_bracket = True
                low_z, low_i, low_f, high_z, high_i, high_f = z, integral, value, z_end, integral_end, value_end
                break
            if z_end >= limit:
                break
            z, decay, temperature, inverse_height, integral, value = z_end, decay_end, t_end, inverse_end, integral_end, value_end
        if not has_bracket:
            height = limit         # the profile does not close inside `MAX_EXOBASE_RADII` (or the panels ran out): the air is taken to reach the limit
            if z_end < limit:
                height = z_end
        else:
            height = high_z
            for _ in range(100):
                if abs(high_f) < abs(low_f):
                    base_z, base_i, base_f = high_z, high_i, high_f
                else:
                    base_z, base_i, base_f = low_z, low_i, low_f
                base_t = temperature_k - rise * math.exp(-base_z * inv_length)
                base_inverse = coefficient / (math.pow(radius_m + base_z, two) * base_t)
                slope = 2.0 / (radius_m + base_z) - base_inverse
                if slope != 0.0:
                    x = base_z - base_f / slope
                else:
                    x = 0.5 * (low_z + high_z)
                if not low_z < x < high_z:
                    x = 0.5 * (low_z + high_z)
                half = 0.5 * (x - base_z)
                middle = base_z + half
                offset = inner * half
                z_a, z_b = middle - offset, middle + offset
                t_a = temperature_k - rise * math.exp(-z_a * inv_length)
                t_b = temperature_k - rise * math.exp(-z_b * inv_length)
                t_x = temperature_k - rise * math.exp(-x * inv_length)
                inverse_x = coefficient / (math.pow(radius_m + x, two) * t_x)
                integral_x = base_i + half * ((base_inverse + inverse_x) / 6.0
                                              + (5.0 / 6.0) * (coefficient / (math.pow(radius_m + z_a, two) * t_a) + coefficient / (math.pow(radius_m + z_b, two) * t_b)))
                value_x = log_offset - integral_x - math.log(t_x) - math.log(inverse_x)
                if value_x > 0.0:
                    low_z, low_i, low_f = x, integral_x, value_x
                else:
                    high_z, high_i, high_f = x, integral_x, value_x
                converged = abs(x - height) <= 1e-12 * x or high_z - low_z <= 1e-12 * high_z or value_x == 0.0
                height = x
                if converged:
                    break
    temperature_at = temperature_k - rise * math.exp(-height / length)
    scale = _K * temperature_at / (over_mass_gravity * math.pow(radius_m / (radius_m + height), two))
    return radius_m + height, 1.0 / (cross_section_m2 * scale), scale


