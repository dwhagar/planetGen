#!/usr/bin/env python3
"""
Fits the three Rayleigh scales of a trojan swarm (GEN.1643, stand-ins TJ9 to TJ11 of docs/design/research/gap-filling-plan.md) to the published Trojan fits.

    trojan_swarm.LIBRATION_SIGMA_DEG   TJ9   libration amplitude of a swarm's bodies
    trojan_swarm.ECC_SIGMA             TJ10  eccentricity
    trojan_swarm.INCLINATION_SIGMA_DEG TJ11  inclination to the secondary's plane

The generator draws each of the three from a Rayleigh distribution. The papers fit other shapes, so the fit is on the shape: a published distribution is turned into its **Rayleigh equivalent**,
the Rayleigh scale that has the same mean (`mean / sqrt(pi/2)`) and, as a check, the same median (`median / sqrt(2 ln 2)`). The two agree when the shape is close to Rayleigh, and the gap
between them says how far it is.

* The Jupiter Trojans ([R5239], the debiased background population of Table 3 for the L4 and L5 clouds) are given as three independent laws in the proper amplitudes A (a sine of inclination),
  B (an eccentricity) and C (the integral of the planar libration, Eq. 3 and 4 of the paper):

      dN/dA ~ A^a1 exp(-1/2 (A/sA)^a2)    dN/dB ~ B^b1 exp(-1/2 (B/sB)^b2)    dN/dC ~ exp(-1/2 ((C - 1.5)/sC)^g)

  A is turned into an inclination with `asin`; C is turned into the libration amplitude (half the peak-to-peak range of the resonant angle, as the Neptune paper [R5238] gives it) by solving
  `f(sigma) = C` of the paper's Eq. 4 for the two ends of the libration, with C cut at 2.5 (the paper's upper bound for a tadpole).
* The Neptune Trojans [R5238] are given directly as Rayleigh widths (Table 4 and the abstract), so their Rayleigh equivalent is the width itself.

It prints the table the report quotes and the proposed values; it changes no constant. Run it with no arguments.

    python3 scripts/fit_trojan_scales.py
"""
import math

C_MIN = 1.5
"""The least value of the libration integral C: the orbit that sits at the point ([R5239: Eq. 3])."""

C_MAX = 2.5
"""The greatest value of C for a tadpole ([R5239: Eq. 3 to 4])."""

A_MAX = 0.7
"""The upper end of the paper's definition interval for A, the proper sine of inclination ([R5239: Table 2])."""

B_MAX = 0.3
"""The upper end taken for B, the proper eccentricity (the paper's stable Trojans lie well inside it)."""

# The debiased background law of each cloud, Table 3 of [R5239]: medians of the posterior. N is the complete population to H 14.5, used as the weight of a cloud.
JUPITER = {
    "L4": {"sA": 0.0533, "a1": 1.274, "a2": 1.008, "sB": 0.0191, "b1": 1.156, "b2": 1.144, "sC": 0.0507, "g": 0.941, "N": 3951},
    "L5": {"sA": 0.3352, "a1": 0.644, "a2": 3.012, "sB": 0.0307, "b1": 1.106, "b2": 1.140, "sC": 0.0443, "g": 0.929, "N": 2664},
}

# The fits to the observed (biased) Trojans, caption of Figure 10 of [R5239]: a second reading of the same shape, kept to show how far the debiasing moves it.
JUPITER_BIASED = {
    "L4": {"sA": 0.095, "a1": 1.0, "a2": 1.3, "sB": 0.04, "b1": 0.87, "b2": 1.75, "sC": 0.045, "g": 1.0},
    "L5": {"sA": 0.15, "a1": 1.2, "a2": 1.6, "sB": 0.038, "b1": 0.85, "b2": 1.6, "sC": 0.045, "g": 1.0},
}

# The Neptune Trojans, [R5238]: the Rayleigh width of the libration amplitude (abstract; Table 4 uses 12), of the eccentricity, and of sin(i) as a width in degrees.
NEPTUNE = {
    "single": {"libration": 15.0, "ecc": 0.045, "inc": 14.0},
    "libration_table4": {"libration": 12.0},
    "cold": {"ecc": 0.02, "inc": 6.0},
    "hot": {"ecc": 0.05, "inc": 18.0},
    "ranges": {"libration": (5.0, 15.0), "ecc": (0.045, 0.045), "inc": (12.0, 16.0)},
}

RAYLEIGH_MEAN = math.sqrt(math.pi / 2.0)
RAYLEIGH_MEDIAN = math.sqrt(2.0 * math.log(2.0))


def generalised(power, scale, shape):
    """The unnormalised density `x**power * exp(-0.5 * (x / scale)**shape)` of the paper's Eq. 6 and 7 (power 0 for Eq. 8)."""
    return lambda x: (x ** power) * math.exp(-0.5 * (x / scale) ** shape)


def moments(pdf, low, high, steps=20000):
    """
    The mean and the median of an unnormalised density on `low`..`high`, by the trapezoid rule on `steps` cells.

    Returns:
        tuple: (mean, median).
    """
    h = (high - low) / steps
    xs = [low + h * (i + 0.5) for i in range(steps)]
    ws = [pdf(x) for x in xs]
    total = sum(ws)
    mean = sum(x * w for x, w in zip(xs, ws)) / total
    run, median = 0.0, high
    for x, w in zip(xs, ws):
        run += w
        if run >= 0.5 * total:
            median = x
            break
    return mean, median


def libration_function(sigma_deg):
    """`f(sigma)` of Eq. 4 of [R5239]: `C = f` on the turning points of the libration of the angle `sigma` from Jupiter."""
    s = abs(math.sin(math.radians(sigma_deg) / 2.0))
    return (1.0 + 4.0 * s ** 3) / (2.0 * s)


def libration_amplitude_deg(c):
    """
    The libration amplitude (half the range of `sigma` between its two turning points) of an orbit of integral `c`, degrees.

    Solves `libration_function(sigma) = c` by bisection on each side of the point at 60 degrees: the near side runs from 60 down to the closest approach (23.9 degrees at C 2.5), the far side from 60
    up to 180 (C 2.5). Zero at `C_MIN`.
    """
    if c <= C_MIN:
        return 0.0
    c = min(c, C_MAX)

    def solve(lo, hi, increasing):
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if (libration_function(mid) < c) == increasing:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    near = solve(1.0, 60.0, False)
    far = solve(60.0, 180.0, True)
    return 0.5 * (far - near)


def separatrix_amplitude_deg():
    """The libration amplitude at the edge of the tadpole zone, `C_MAX`: half the range of `sigma` between the closest approach to the secondary and the far side (78 degrees)."""
    return libration_amplitude_deg(C_MAX)


def closest_approach_deg():
    """The smallest `sigma` a tadpole reaches (the near turning point at `C_MAX`): the 23.9 degrees that `trojans.TADPOLE_LIMIT_DEG` holds, which is a distance from the secondary, not an amplitude."""
    lo, hi = 1.0, 60.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if libration_function(mid) > C_MAX:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def truncated_rayleigh_mean(sigma, limit):
    """The mean of a Rayleigh variate of scale `sigma` cut at `limit` (the closed form of what `trojan_swarm.truncated_rayleigh` draws)."""
    top = 1.0 - math.exp(-0.5 * (limit / sigma) ** 2)
    # integral of x * x/s^2 exp(-x^2/2s^2) from 0 to limit = s sqrt(pi/2) erf(limit / (s sqrt 2)) - limit exp(-limit^2 / 2 s^2)
    mean = sigma * RAYLEIGH_MEAN * math.erf(limit / (sigma * math.sqrt(2.0))) - limit * math.exp(-0.5 * (limit / sigma) ** 2)
    return mean / top


def share_beyond(params, amplitude_deg):
    """The share of a cloud's C law whose libration amplitude exceeds `amplitude_deg`."""
    pdf = generalised(0.0, params["sC"], params["g"])
    steps = 4000
    h = (C_MAX - C_MIN) / steps
    cs = [C_MIN + h * (i + 0.5) for i in range(steps)]
    ws = [pdf(c - C_MIN) for c in cs]
    total = sum(ws)
    return sum(w for c, w in zip(cs, ws) if libration_amplitude_deg(c) > amplitude_deg) / total


def amplitude_pair(params):
    """The mean and median of the libration amplitude (degrees) of a cloud's C law, cut at `C_MAX`."""
    pdf = generalised(0.0, params["sC"], params["g"])
    steps = 4000
    h = (C_MAX - C_MIN) / steps
    cs = [C_MIN + h * (i + 0.5) for i in range(steps)]
    ws = [pdf(c - C_MIN) for c in cs]
    amps = [libration_amplitude_deg(c) for c in cs]
    total = sum(ws)
    mean = sum(a * w for a, w in zip(amps, ws)) / total
    run, median = 0.0, amps[-1]
    for a, w in zip(amps, ws):
        run += w
        if run >= 0.5 * total:
            median = a
            break
    return mean, median


def inclination_pair(params):
    """The mean and median of the inclination (degrees, `asin` of A) of a cloud's A law."""
    mean_a, median_a = moments(generalised(params["a1"], params["sA"], params["a2"]), 0.0, A_MAX)
    # the mean of asin(A) is not asin of the mean: integrate the angle directly
    pdf = generalised(params["a1"], params["sA"], params["a2"])
    steps = 20000
    h = A_MAX / steps
    xs = [h * (i + 0.5) for i in range(steps)]
    ws = [pdf(x) for x in xs]
    total = sum(ws)
    mean = sum(math.degrees(math.asin(x)) * w for x, w in zip(xs, ws)) / total
    return mean, math.degrees(math.asin(median_a))


def eccentricity_pair(params):
    """The mean and median of a cloud's B law."""
    return moments(generalised(params["b1"], params["sB"], params["b2"]), 0.0, B_MAX)


def equivalents(pair):
    """The Rayleigh equivalent by the mean and by the median of a (mean, median) pair."""
    return pair[0] / RAYLEIGH_MEAN, pair[1] / RAYLEIGH_MEDIAN


def fit(table):
    """
    The Rayleigh equivalents of one set of cloud laws.

    Returns:
        dict: For each cloud and for the N-weighted pair of clouds (when `N` is given), `{"libration": (mean_eq, median_eq), "ecc": ..., "inc": ...}` plus the raw (mean, median).
    """
    out = {}
    for cloud, params in table.items():
        pairs = {"libration": amplitude_pair(params), "ecc": eccentricity_pair(params), "inc": inclination_pair(params)}
        out[cloud] = {key: {"mean": p[0], "median": p[1], "eq_mean": equivalents(p)[0], "eq_median": equivalents(p)[1]} for key, p in pairs.items()}
    weights = {cloud: params.get("N") for cloud, params in table.items()}
    if all(weights.values()):
        total = sum(weights.values())
        both = {}
        for key in ("libration", "ecc", "inc"):
            both[key] = {field: sum(weights[c] * out[c][key][field] for c in table) / total for field in ("mean", "median", "eq_mean", "eq_median")}
        out["L4+L5"] = both
    return out


def report():
    """The lines the report quotes."""
    lines = []
    for name, table in (("Jupiter, debiased background (Table 3)", JUPITER), ("Jupiter, observed (Figure 10 caption)", JUPITER_BIASED)):
        lines.append(name)
        for cloud, got in fit(table).items():
            for key, unit in (("libration", "deg"), ("ecc", ""), ("inc", "deg")):
                g = got[key]
                lines.append(f"  {cloud:6s} {key:9s} mean {g['mean']:8.4f}  median {g['median']:8.4f}  Rayleigh scale by mean {g['eq_mean']:8.4f}  by median {g['eq_median']:8.4f} {unit}")
    lines.append("The tadpole zone (Eq. 3 and 4 of the Jupiter paper)")
    lines.append(f"  closest approach of a tadpole {closest_approach_deg():.2f} deg from the secondary; amplitude at the edge of the zone {separatrix_amplitude_deg():.2f} deg")
    for cloud, params in JUPITER.items():
        lines.append(f"  {cloud}: share of the Jupiter law with amplitude over 23.9 deg {share_beyond(params, 23.9):.3f}, over 36.1 deg {share_beyond(params, 36.1):.3f}")
    for sigma in (10.0, 13.0):
        lines.append(f"  mean of the generator's truncated draw at scale {sigma:.0f}: cut 23.9 deg {truncated_rayleigh_mean(sigma, 23.9):.2f}, cut {separatrix_amplitude_deg():.1f} deg "
                     f"{truncated_rayleigh_mean(sigma, separatrix_amplitude_deg()):.2f}")
    lines.append("Neptune, Rayleigh widths (Table 4 and abstract)")
    for label, row in NEPTUNE.items():
        lines.append(f"  {label:17s} {row}")
    return lines


def main():
    print("\n".join(report()))


if __name__ == "__main__":
    main()
