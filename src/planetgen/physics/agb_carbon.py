# planetgen/physics/agb_carbon.py

"""
Carbon stars and S stars from the thermal pulses of the asymptotic giant branch (GEN.1026, after GEN.944).

The post-main-sequence track (`post_main_sequence`) has no thermal pulses: its asymptotic branch is one stretch. This module counts them from the
star's place in that stretch and lets the third dredge-up carry carbon to the surface, so a star's surface C/O follows from its mass, its [Fe/H]
and how far along the branch it is. It is a derived value, as `variability.kind_of` is: the sheet computes it from the stored initial mass, mass,
age, [Fe/H] and `state`, and nothing is stored.

**Pulses.** `pulses_in_branch` is the number of thermal pulses over the whole branch for a star of the initial mass: the table of [R2244] (Table 1, Z 0.007 and 0.014, 17 to 25 pulses from 1.5 to 3.75 solar masses). The star has run that many times its share of the
branch's age window.

**Dredge-up.** From the first pulse a star above `tdu_minimum_mass_sol` brings carbon up: at pulse `n` the dredged share of the core growth is
`lambda(n) = lambda_max (1 - exp(-n / DREDGE_RISE_PULSES))` (the form of [R2243] as written in the population-synthesis papers),
`lambda_max` the table of [R2244] by mass. Each pulse grows the core by `CORE_GROWTH_PER_PULSE_SOL`; the dredged mass is
`lambda` times that, of which `INTERSHELL_CARBON_FRACTION` is carbon. The surface C/O is that carbon (and the envelope's own) over the envelope's
oxygen, with the envelope the mean of its mass at the start of the branch and its mass now, less the core. Above `HBB_MIN_MASS_SOL` the base of the
envelope burns the carbon again (hot bottom burning) and the star stays oxygen-rich.

A star on the asymptotic branch with C/O at or above 1 is a **carbon star**; one with C/O from `S_STAR_MIN_C_TO_O` to 1 whose dredge-up has run is an
**S star** (the zirconium oxide bands fade as the oxygen is taken up by carbon; the s-process products come up with the carbon).
"""

import math

from planetgen import tuning
from planetgen.physics import post_main_sequence

CARBON_STAR = "carbon star"
S_STAR = "S star"
KINDS = (CARBON_STAR, S_STAR)

PULSES_BY_MASS_SOL = ((1.5, 17.0), (1.75, 19.0), (2.75, 25.0), (3.25, 22.0), (3.75, 22.5))
"""tuple: `(initial mass, thermal pulses over the asymptotic branch)`, linear between, held outside. [Stand-in CS1, tier T3, tweak L0: [R2244: Table 1] ([R2244]; Z 0.007 and 0.014 rows; 1.5 solar masses has 18 and 16, 1.75 has 17 and 20, 2.75 25, 3.25 22, 3.75 22 and 23).]"""

LAMBDA_MAX_BY_MASS_SOL = ((1.5, 0.41), (1.75, 0.50), (2.75, 0.80), (3.25, 0.92), (3.75, 0.96))
"""tuple: `(initial mass, largest dredge-up share of the core growth)`, linear between, held outside. [Stand-in CS2, tier T3, tweak L0: the
lambda_max column of [R2244: Table 1] (1.5: 0.31 at Z 0.007 and 0.51 at 0.014; 1.75: 0.51 and 0.50; 2.75: 0.80; 3.25: 0.92; 3.75: 0.97
and 0.95).]"""

DREDGE_RISE_PULSES = 4.0
"""float: Pulses over which the dredge-up share rises to 63 percent of its largest. [Stand-in CS3, tier T6, tweak L1: [R2243] fit it to
models; their coefficients are not at hand, and a handful of pulses is what the published lambda-against-pulse curves show.]"""

CORE_GROWTH_PER_PULSE_SOL = 0.004
"""float: Core mass added by one thermal pulse cycle, solar masses. [Stand-in CS4, tier T6, tweak L2: the core grows by about 0.1 solar masses over
the 20 to 25 pulses of the models in CS1, 0.004 per pulse.]"""

INTERSHELL_CARBON_FRACTION = 0.22
"""float: The carbon mass fraction of the dredged intershell material. [Stand-in CS5, tier T6, tweak L1: about a fifth to a quarter of
the intershell is carbon in the published AGB models; read second-hand.]"""

CORE_MASS_SOL = 0.6
"""float: The core's mass for the envelope, solar masses. [Stand-in CS6, tier T6, tweak L0: the helium core of a thermally pulsing star is 0.55
to 0.7.]"""

TDU_MIN_MASS_AT_SOLAR_SOL = 1.5
"""float: The lightest star at solar [Fe/H] that dredges carbon up. [Stand-in CS7, tier T3, tweak L1: [R2244] (1.5 solar masses
becomes carbon-rich at Z 0.007 and 0.014 in their models with overshoot); [R0925] find 1.5 at Z 0.02 and 1.4 at 0.01; the models
without overshoot put the solar minimum at 2.]"""

TDU_MIN_MASS_SLOPE_PER_DEX = 0.2
"""float: How the lightest dredging mass falls per dex of [Fe/H] below solar, solar masses (a metal-poor star dredges at lower mass). [Stand-in
CS8, tier T6, tweak L2: the Cristallo step from 1.5 at Z 0.02 to 1.4 at 0.01 is 0.33 per dex and the carbon stars of the metal-poor dwarf galaxies
reach below 1.2, so 0.2 per dex.]"""

TDU_MIN_MASS_FLOOR_SOL = 1.0
"""float: The lightest star that dredges at any [Fe/H]. [Stand-in CS9, tier T6, tweak L1.]"""

HBB_MIN_MASS_SOL = 4.0
"""float: A star of this mass or more burns its dredged carbon at the base of its envelope and stays oxygen-rich. [Stand-in CS10, tier T3, tweak L0:
[R2244] (hot bottom burning from about 4 solar masses); Izzard-style synthesis [R0796] takes the same 4.0 as the
upper edge of the carbon stars.]"""

S_STAR_MIN_C_TO_O = 0.9
"""float: A star with C/O from this to 1 whose dredge-up has run is an S star. [Stand-in CS11, tier T6, tweak L1: S stars run from about 0.5 to
0.99 in C/O but the ones with the dredged s-process come up with the carbon, so near 1.]"""


def _interpolated(table, mass_sol):
    if mass_sol <= table[0][0]:
        return table[0][1]
    for (low, a), (high, b) in zip(table, table[1:]):
        if mass_sol <= high:
            return a + (b - a) * (mass_sol - low) / (high - low)
    return table[-1][1]


def pulses_in_branch(initial_mass_sol):
    """Thermal pulses over the whole asymptotic branch of a star of this initial mass (CS1)."""
    return _interpolated(PULSES_BY_MASS_SOL, initial_mass_sol)


def tdu_minimum_mass_sol(fe_h):
    """The lightest star of birth [Fe/H] `fe_h` that dredges carbon up (CS7 to CS9)."""
    return max(TDU_MIN_MASS_AT_SOLAR_SOL + TDU_MIN_MASS_SLOPE_PER_DEX * min(fe_h, 0.0), TDU_MIN_MASS_FLOOR_SOL)


def dredged_core_mass_sol(initial_mass_sol, pulses):
    """The core growth brought up over `pulses` thermal pulses (a whole number), solar masses: the sum over the pulses of `lambda(n)` times
    the core growth of one (CS2 to CS4)."""
    lambda_max = _interpolated(LAMBDA_MAX_BY_MASS_SOL, initial_mass_sol)
    total = sum(lambda_max * (1.0 - math.exp(-n / DREDGE_RISE_PULSES)) for n in range(1, pulses + 1))
    return total * CORE_GROWTH_PER_PULSE_SOL


def _solar_mass_fraction(element, atomic_mass):
    """The envelope mass fraction of `element` at solar [Fe/H]: A21's log eps against the proto-solar hydrogen (`tuning`)."""
    return 10.0 ** (tuning.STELLAR_SOLAR_LOG_EPS[element] - 12.0) * atomic_mass * tuning.STELLAR_PROTO_X / 1.008


def surface_c_to_o(initial_mass_sol, fe_h, mass_start_sol, mass_now_sol, pulses):
    """
    The surface carbon-to-oxygen number ratio of a star of birth [Fe/H] `fe_h` after `pulses` thermal pulses, with the envelope the mean of the
    mass it had when the branch began (`mass_start_sol`) and has now (`mass_now_sol`) less the core; the ratio before any dredge-up (its own, from
    the A21 abundances) when the star is too light to dredge or has run no pulse.
    """
    scale = 10.0 ** fe_h
    carbon = _solar_mass_fraction("C", 12.011) * scale / 12.011
    oxygen = _solar_mass_fraction("O", 15.999) * scale / 15.999
    if pulses < 1 or initial_mass_sol < tdu_minimum_mass_sol(fe_h):
        return carbon / oxygen
    envelope = max(0.5 * (mass_start_sol + mass_now_sol) - CORE_MASS_SOL, 0.05)
    added = INTERSHELL_CARBON_FRACTION * dredged_core_mass_sol(initial_mass_sol, pulses) / 12.011 / envelope
    return (carbon + added) / oxygen


def initial_mass_from_lifespan_sol(lifespan_gy, fe_h):
    """The initial mass, solar masses, of the star of birth [Fe/H] `fe_h` that becomes a remnant at `lifespan_gy` (the end of its last giant
    phase, the `lifespan_gy` a star stores): the track's end age falls with mass, so a bisection finds it. None when no star of 0.8 to 8 solar
    masses has that lifespan (a lifespan that is not a star's)."""
    low, high = 0.8, 8.0
    if lifespan_gy is None or not post_main_sequence.timeline(high, fe_h).end_gy <= lifespan_gy <= post_main_sequence.timeline(low, fe_h).end_gy:
        return None
    for _ in range(40):
        middle = math.sqrt(low * high)
        if post_main_sequence.timeline(middle, fe_h).end_gy > lifespan_gy:
            low = middle
        else:
            high = middle
    return math.sqrt(low * high)


def pulses_run(initial_mass_sol, fe_h, age_gy):
    """The thermal pulses a star of this initial mass and birth [Fe/H] has run by `age_gy`: its share of the asymptotic window (the track's) times
    the pulses of the whole branch, rounded down; 0 outside the window."""
    window = post_main_sequence.timeline(initial_mass_sol, fe_h).window(post_main_sequence.ASYMPTOTIC_GIANT)
    if window is None or not window.start_gy <= age_gy < window.end_gy:
        return 0
    share = (age_gy - window.start_gy) / (window.end_gy - window.start_gy)
    return int(pulses_in_branch(initial_mass_sol) * share)


def kind_of(state, initial_mass_sol, mass_sol, fe_h, age_gy):
    """
    `CARBON_STAR`, `S_STAR` or None for a star of this `state` (`stellar_state.STATES`), initial and present mass (solar masses), birth [Fe/H]
    and age (Gy): only a star on the asymptotic giant branch, lighter than `HBB_MIN_MASS_SOL`, that has dredged carbon up is either.
    """
    if state != "asymptotic-giant" or None in (initial_mass_sol, mass_sol, fe_h, age_gy) or initial_mass_sol >= HBB_MIN_MASS_SOL:
        return None
    if initial_mass_sol < tdu_minimum_mass_sol(fe_h):
        return None
    pulses = pulses_run(initial_mass_sol, fe_h, age_gy)
    if pulses < 1:
        return None
    window = post_main_sequence.timeline(initial_mass_sol, fe_h).window(post_main_sequence.ASYMPTOTIC_GIANT)
    mass_start = post_main_sequence.timeline(initial_mass_sol, fe_h).state(window.start_gy).mass_sol
    ratio = surface_c_to_o(initial_mass_sol, fe_h, mass_start, mass_sol, pulses)
    if ratio >= 1.0:
        return CARBON_STAR
    if ratio >= S_STAR_MIN_C_TO_O:
        return S_STAR
    return None
