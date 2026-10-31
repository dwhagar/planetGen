# planetgen/physics/impact_spin.py

"""
The birth spin of a rocky planet from the mass ratio of its last impacts (GEN.978)
==================================================================================

Before GEN.979 `spin.accretion_spin` drew a rocky planet's birth spin as three Gaussian components of one width fixed by Kokubo and Genda 2010 (ApJ 714:L21,
arXiv:1003.4384): 0.41 of the breakup spin per component, whatever the planet's mass, which put every planet at a 2 to 3 hour day. That width is the
average over the planets of their N-body runs, whose last impacts were collisions between protoplanets of comparable mass. The paper's equation 4 says what
sets the spin: the angular momentum the last one or few collisions leave,

    omega ~ (5 / sqrt 2) <g(gamma)^2>^(1/2) <sin^2 theta>^(1/2) <nu^2>^(1/2) omega_cr,
    g(gamma) = gamma (1 - gamma) [gamma^(1/3) + (1 - gamma)^(1/3)]^(1/2),

with `gamma = m_impactor / (m_impactor + m_target)`, `theta` the angle of the impact, `nu = v_impact / v_escape` and `omega_cr = sqrt(G M / R^3)` the breakup spin
(read in full text, 10-29). A planet hit by a body of its own size has `g = 0.315` and spins near breakup; one built from small bodies (`g` is about `gamma`)
spins slowly, which is how Mars (24.6 h, 0.07 of breakup) comes about without a fitted scale. This module is that law and a **stand-in for the
distribution of `gamma` by planet mass**, which the papers do not publish in a form we hold (Kokubo and Ida 2007 and Agnor, Canup and Levison 1999 hold
the N-body statistics; Boss accepted a registered stand-in for it, 10-29 8:48 PM EDT). `spin.accretion_spin` calls it (GEN.979).

**The vector rms.** With the impact parameter uniform in area, `<sin^2 theta> = 1/2`; with the impact velocity a little above the escape velocity,
`nu_rms = IMPACT_SPEED_RATIO` (1.1). Then `omega_rms / omega_cr = C g(gamma)`, `C = (5/sqrt 2) sqrt(1/2) nu_rms = 2.75`, and each of the three Gaussian
components has `1/sqrt 3` of it. At the `gamma` of 0.28 this is the paper's realistic-model rms of 0.71 breakup (2.35 rad/h at 3.3 rad/h), so
the paper's own runs are protoplanet-on-protoplanet impacts of mass ratio near 0.3, a check this module's tests make.

**The distribution of gamma (stand-in).** A planet either had a giant impact of embryos near its end or was built by planetesimals alone:

* The chance of a giant impact is `1 / (1 + (M_iso / M)^2)`, `M` the planet's mass and `M_iso = ISOLATION_MASS_EARTH` (0.1 Earth masses, the isolation
  mass of oligarchic growth at 1 AU, Kokubo and Ida 2000): well under it a body is an embryo that never met an equal, well over it the planet is made of them.
* With a giant impact `gamma` is log-uniform over `GIANT_GAMMA_RANGE` (0.02 to 0.5; the last impactors of Agnor et al. 1999 and Kokubo and Ida 2007 run from
  a few percent of the planet to its own size).
* Without one it is log-uniform over `PLANETESIMAL_GAMMA_RANGE` (0.001 to 0.03: the effective ratio of the planetesimals that spun the embryo up).

The resulting periods run from 1.3 hours (a near-breakup spin) to some thousand hours: the widest range of the two sources in question, Kokubo and Genda's 2 to 3 hours and Miguel and Brunini
2010's 10 to 10,000 hours (arXiv:1004.1406, read in full), as Boss's rule for two sources in question asks. The tilt is `spin.accretion_spin`'s: isotropic after a giant impact (the paper), prograde on the planetesimal branch (Schlichting and Sari 2007). The stand-ins are in `docs/design/research/gap-filling-plan.md` section 6.26.
"""

import math

IMPACT_SPEED_RATIO = 1.1
"""float: `nu_rms = <(v_impact / v_escape)^2>^(1/2)` of the last impacts: a body falling in from rest at infinity arrives at the escape speed (1), the
protoplanets of an oligarchy a little faster. Source tier T6 (physics fixes the floor at 1, the excess is argued), tweak level L2 (the value that gives the
paper's own rms of 0.71 breakup at the paper's `gamma` of 0.28)."""
MEAN_SIN2_THETA = 0.5
"""float: `<sin^2 theta>` for an impact parameter uniform in area (`b = sin theta`, `p(b) = 2 b`). Source tier T5, tweak level L0."""
ISOLATION_MASS_EARTH = 0.1
"""float: The mass below which a body is an embryo that never had a giant impact, Earth masses (the isolation mass of oligarchic growth at 1 AU, 0.1
to 0.2 in Kokubo and Ida 2000, recalled; the planets of the same disc past it are made by giant impacts). Source tier T6, tweak level L3 (chosen so
that a Mars-sized planet is as likely to have had one as not)."""
ISOLATION_STEEPNESS = 2.0
"""float: The exponent of the logistic in `log M` that gives the chance of a giant impact. Source tier T6, tweak level L3."""
GIANT_GAMMA_RANGE = (0.02, 0.5)
"""tuple: The mass ratio of the last giant impact, log-uniform (impactors of a few percent of the planet to its own size; Agnor et al. 1999, Kokubo and Ida 2007,
recalled; `g` is symmetric about 1/2, so the larger of the two bodies is the target). Source tier T6, tweak level L3 (the lower end is chosen so that the median free
day of 0.1 to 5 Earth masses is 7 hours, the Vesta, Ceres and Earth-Moon range)."""
PLANETESIMAL_GAMMA_RANGE = (0.001, 0.03)
"""tuple: The effective `gamma` of an embryo built by planetesimals alone, log-uniform (spins from a tenth of a percent to 3 percent of breakup per unit
of `C`; Mars at 0.07 of breakup needs 0.03). Source tier T6, tweak level L3 (the upper end is Mars; the lower end gives the slow spins of Miguel and Brunini 2010)."""


def g_function(gamma):
    """`g(gamma)` of the paper's equation 4, symmetric under `gamma -> 1 - gamma`, 0.315 at 1/2."""
    if not 0.0 <= gamma <= 1.0:
        raise ValueError(f"a mass ratio of {gamma} is not between 0 and 1")
    return gamma * (1.0 - gamma) * math.sqrt(gamma ** (1.0 / 3.0) + (1.0 - gamma) ** (1.0 / 3.0))


def spin_constant():
    """`C = (5 / sqrt 2) <sin^2 theta>^(1/2) nu_rms`, the factor of `omega_rms / omega_cr = C g(gamma)`."""
    return 5.0 / math.sqrt(2.0) * math.sqrt(MEAN_SIN2_THETA) * IMPACT_SPEED_RATIO


def rms_breakup_fraction(gamma):
    """The rms of the birth-spin vector in breakup spins for last impacts of mass ratio `gamma`."""
    return spin_constant() * g_function(gamma)


def component_breakup_fraction(gamma):
    """The Gaussian width of each of the vector's three components, in breakup spins: the rms over `sqrt 3`."""
    return rms_breakup_fraction(gamma) / math.sqrt(3.0)


def giant_impact_probability(mass_earth):
    """The chance that a planet of `mass_earth` Earth masses had a giant impact of embryos at its end: `1 / (1 + (M_iso / M)^k)`."""
    if mass_earth <= 0.0:
        raise ValueError(f"a mass of {mass_earth} Earth masses is not a planet")
    return 1.0 / (1.0 + (ISOLATION_MASS_EARTH / mass_earth) ** ISOLATION_STEEPNESS)


def had_giant_impact(mass_earth, branch_u):
    """Whether a planet of `mass_earth` Earth masses is on the giant-impact branch for the uniform number `branch_u`."""
    return branch_u < giant_impact_probability(mass_earth)


def _log_uniform(bounds, u):
    low, high = bounds
    return low * (high / low) ** u


def sample_gamma(mass_earth, branch_u, within_u):
    """
    `gamma` for a planet of `mass_earth` Earth masses from two uniform numbers in `[0, 1)`: `branch_u` below `giant_impact_probability` picks the giant
    branch, and `within_u` places it in that branch's log-uniform range. The caller supplies the numbers (so the random stream is its to place).
    """
    if had_giant_impact(mass_earth, branch_u):
        return _log_uniform(GIANT_GAMMA_RANGE, within_u)
    return _log_uniform(PLANETESIMAL_GAMMA_RANGE, within_u)


def _mean_g2(bounds, steps=400):
    """`<g^2>` over a log-uniform range of gamma, by the midpoint rule in `ln gamma`."""
    low, high = (math.log(b) for b in bounds)
    return sum(g_function(math.exp(low + (high - low) * (i + 0.5) / steps)) ** 2 for i in range(steps)) / steps


_MEAN_G2_GIANT = _mean_g2(GIANT_GAMMA_RANGE)
_MEAN_G2_PLANETESIMAL = _mean_g2(PLANETESIMAL_GAMMA_RANGE)


def mean_rms_breakup_fraction(mass_earth):
    """
    The rms birth spin of planets of `mass_earth` Earth masses averaged over the distribution of `gamma`, in breakup spins: `C sqrt(<g^2>)` with
    the mean over both branches. The width of a draw that does not choose `gamma` for each planet.
    """
    p = giant_impact_probability(mass_earth)
    return spin_constant() * math.sqrt(p * _MEAN_G2_GIANT + (1.0 - p) * _MEAN_G2_PLANETESIMAL)


def component_sigma_rad_s(critical_spin_rad_s, gamma):
    """The Gaussian width of each component of a planet's birth spin vector, rad/s, for last impacts of mass ratio `gamma`: the share
    `component_breakup_fraction` of the body's breakup spin `critical_spin_rad_s`."""
    return component_breakup_fraction(gamma) * critical_spin_rad_s
