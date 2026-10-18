# planetgen/physics/wolf_rayet.py

"""
The Wolf-Rayet state of a very massive star (GEN.243; docs/design/intergalactic-and-halo-contents.md 5.11).

A star of at least `threshold_mass_sol` initial mass loses its hydrogen envelope to its wind and ends its supergiant
phase as a Wolf-Rayet star: a hot helium (WN, nitrogen-rich), then carbon (WC) and, for the heaviest, oxygen (WO)
core with a dense wind, for `tuning.WR_PHASE_MYR` before it collapses. The state replaces the end of the supergiant
window of `stellar_evolution.evolve_star`, so the bright and mass passes that place the massive stars draw Wolf-Rayet
stars with no rate of their own.

Inputs that are inferred, not sourced, are flagged [I] in `tuning`: the oxygen-rich threshold, the temperature split
by subtype, and the exponent of the metallicity dependence, which is set so the thresholds at half and a fifth of the
Sun's metallicity come out at 35 and 56 solar masses.
"""

from planetgen import tuning

CLASS = "W"
"""str: The `yerkes_class` of a Wolf-Rayet star."""
CLASS_NAME = "Wolf-Rayet"
"""str: Display name of the class, in the place `YERKES_CLASS_NAMES` is used."""
MYR_PER_GY = 1.0e3


def threshold_mass_sol(fe_h=None):
    """The least initial mass that ends as a Wolf-Rayet star at this [Fe/H] (the Sun's when None), solar masses."""
    metals = 10.0 ** (fe_h or 0.0)
    return tuning.WR_MIN_INITIAL_MASS_SOL * metals ** -tuning.WR_THRESHOLD_METALLICITY_EXPONENT


def phase_gy(supergiant_window_gy):
    """How long the phase lasts in a star whose supergiant window is this long, Gy."""
    return min(tuning.WR_PHASE_MYR / MYR_PER_GY, tuning.WR_MAX_SUPERGIANT_SHARE * supergiant_window_gy)


def start_gy(mass_sol, main_sequence_gy, end_gy, fe_h=None):
    """
    The age a star of this initial mass becomes Wolf-Rayet, Gy, or None when it never does. `end_gy` is where its
    supergiant window (and so its life) ends.
    """
    if mass_sol < threshold_mass_sol(fe_h):
        return None
    return end_gy - phase_gy(end_gy - main_sequence_gy)


def subtype(mass_sol, progress):
    """`"WN"`, `"WC"` or `"WO"` for a star of this initial mass `progress` (0 to 1) of the way through its phase."""
    gone = 0.0
    for name, share in tuning.WR_SUBTYPE_PHASE_SHARES:
        gone += share
        if progress < gone:
            break
    if name == "WO" and mass_sol < tuning.WO_MIN_INITIAL_MASS_SOL:
        return "WC"
    return name


def subclass(name, temperature_k):
    """The subclass number of a star of this subtype and effective temperature: the hotter, the lower, over the
    subtype's range of both."""
    low_k, high_k = tuning.WR_TEMPERATURE_RANGE_K[name]
    cool, hot = tuning.WR_SUBCLASS_RANGE[name][1], tuning.WR_SUBCLASS_RANGE[name][0]
    share = min(max((temperature_k - low_k) / (high_k - low_k), 0.0), 1.0)
    return round(cool + (hot - cool) * share)


def companion_chance(star_type):
    """The share of stars with this spectral text (`"WN6 Blue Wolf-Rayet Star"`) that have a companion
    (`tuning.WR_COMPANION_CHANCE`), or None for a star that is not a Wolf-Rayet star, which keeps the field's rate."""
    return tuning.WR_COMPANION_CHANCE.get((star_type or "")[:2])


def drawn_type(rng=None):
    """The spectral text of a Wolf-Rayet star drawn with nothing else known (`"WN6"`): a subtype by the phase shares,
    a heavy star's, and a subclass uniformly in its range."""
    from planetgen.util import draw
    rng = rng or draw
    name = subtype(tuning.WO_MIN_INITIAL_MASS_SOL, rng.random())
    low, high = sorted(tuning.WR_SUBCLASS_RANGE[name])
    return f"{name}{rng.randint(low, high)}"


def state(mass_sol, age_gy, start, end, rng, min_luminosity_sol=None):
    """
    The state of a Wolf-Rayet star of initial mass `mass_sol` at `age_gy`, between `start` and `end` (Gy).

    Returns:
        dict: `yerkes_class` (`CLASS`), `mass_sol` (the present mass), `luminosity_sol`, `temperature_k`,
            `wr_subtype` (`"WN"`, `"WC"` or `"WO"`) and `phase_end_gy`; a luminosity floor `min_luminosity_sol`
            (the bright pass's) draws from the part of the range at or above it.
    """
    progress = (age_gy - start) / (end - start)
    name = subtype(mass_sol, progress)
    low, high = tuning.WR_LUMINOSITY_RANGE_SOL
    if min_luminosity_sol is not None:
        low = min(max(low, min_luminosity_sol), high)
    luminosity = low * (high / low) ** rng.random()
    return {"yerkes_class": CLASS, "wr_subtype": name,
            "mass_sol": mass_sol * rng.uniform(*tuning.WR_PRESENT_MASS_FRACTION_RANGE),
            "luminosity_sol": luminosity, "temperature_k": rng.uniform(*tuning.WR_TEMPERATURE_RANGE_K[name]),
            "phase_end_gy": end}
