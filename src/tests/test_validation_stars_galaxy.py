# tests/test_validation_stars_galaxy.py

"""GEN.753: the statistics of the star function and the stellar-to-halo relation that agree with their published sources
(project report gen753-stars-galaxy-validation-2026-10-31.md). Analytic, no dataset and no draws: each test fails if a later change
moves the generator outside its source's range. The Milky Way's star formation rate is a test target here, never an input."""

import math

from planetgen import tuning
from planetgen.galaxy import mass_chain
from planetgen.physics import imf

STAR_FORMATION_SOL_PER_YR = 1.9
"""float: The Galaxy's star formation rate, Chomiuk and Povich 2011 (quoted in Rosslowe and Crowther 2015 section 4.4)."""
H = 0.671
OMEGA_MATTER = 0.3175


def _per_solar_mass_formed(low_sol):
    return imf.number(low_sol, imf.highest_mass_sol()) / imf.integral(1.0, imf.lowest_mass_sol(), imf.highest_mass_sol())


def _chabrier(mass_sol):
    """Chabrier 2003, individual stars of the disk, per unit of ln mass (any scale; parameters recalled)."""
    inside = math.exp(-(math.log10(min(mass_sol, 1.0)) - math.log10(0.079)) ** 2 / (2 * 0.69 ** 2))
    return inside * (mass_sol ** -1.3 if mass_sol > 1.0 else 1.0)


def test_imf_shares_match_chabrier_2003_in_seven_mass_bins():
    low, high = imf.lowest_mass_sol(), imf.highest_mass_sol()
    steps = 20000
    edges = [low * (high / low) ** (i / steps) for i in range(steps + 1)]
    mids = [math.sqrt(edges[i] * edges[i + 1]) for i in range(steps)]
    weights = [_chabrier(m) * math.log(edges[i + 1] / edges[i]) for i, m in enumerate(mids)]
    total = sum(weights)
    for a, b in ((0.08, 0.2), (0.2, 0.5), (0.5, 1), (1, 2), (2, 8), (8, 25), (25, 150)):
        expected = sum(w for m, w in zip(mids, weights) if a <= m < b) / total
        got = imf.number(a, b) / imf.number(low, high)
        assert abs(got - expected) <= 0.01 * expected, (a, b, got, expected)


def test_the_galactic_wolf_rayet_count_is_inside_rosslowe_and_crowther_2015():
    """1200 (+300, -100) from star formation 1.9 solar masses a year, progenitors above 25 solar masses (section 4.4)."""
    per_solar_mass = _per_solar_mass_formed(tuning.WR_MIN_INITIAL_MASS_SOL)
    count = per_solar_mass * STAR_FORMATION_SOL_PER_YR * 1.0e6 * tuning.WR_PHASE_MYR
    assert 1100.0 <= count <= 1500.0, count


def test_the_core_collapse_rate_is_inside_diehl_2006_from_26_aluminium():
    """1.9 +/- 1.1 a century (Nature 439, 45)."""
    rate = _per_solar_mass_formed(tuning.SUPERGIANT_MIN_MASS_SOL) * STAR_FORMATION_SOL_PER_YR * 100.0
    assert 0.8 <= rate <= 3.0, rate


# Leauthaud et al. 2012 Eq. 13, Table 5 SIG_MOD1 z1: halos at 200 times the mean density (M200b)
_LEAUTHAUD = dict(log_m1=12.520, log_ms0=10.916, beta=0.457, delta=0.566, gamma=1.53)
# Behroozi et al. 2019 Eq. J1, Table J1 (observed stellar mass, central galaxies, no intrahalo light, z = 0): virial halos
_BEHROOZI = dict(eps=-1.435, m1=12.081, alpha=1.957, beta=0.474, delta=0.386, log_gamma=-1.065)


def _nfw(x):
    return math.log(1.0 + x) - x / (1.0 + x)


def _mass_at_overdensity(m200c, overdensity_critical):
    """The mass inside the radius of this mean density (times critical) of an NFW halo of this M200c (Dutton and Maccio
    2014 concentration at z = 0)."""
    c = 10 ** (0.905 - 0.101 * math.log10(m200c / (1e12 / H)))
    low, high = 0.05, 20.0
    for _ in range(100):
        x = math.sqrt(low * high)
        if 200.0 * _nfw(c * x) / _nfw(c) / x ** 3 > overdensity_critical:
            low = x
        else:
            high = x
    x = math.sqrt(low * high)
    return m200c * _nfw(c * x) / _nfw(c)


def _leauthaud_log_halo(log_stellar):
    p = _LEAUTHAUD
    ratio = 10 ** (log_stellar - p["log_ms0"])
    return p["log_m1"] + p["beta"] * (log_stellar - p["log_ms0"]) + ratio ** p["delta"] / (1 + ratio ** -p["gamma"]) - 0.5


def _leauthaud_log_stellar(log_m200b):
    low, high = 6.0, 13.5
    for _ in range(80):
        middle = 0.5 * (low + high)
        low, high = (middle, high) if _leauthaud_log_halo(middle) < log_m200b else (low, middle)
    return 0.5 * (low + high)


def _behroozi_log_stellar(log_mvir):
    p = _BEHROOZI
    x = log_mvir - p["m1"]
    return p["m1"] + p["eps"] - math.log10(10 ** (-p["alpha"] * x) + 10 ** (-p["beta"] * x)) \
        + 10 ** p["log_gamma"] * math.exp(-0.5 * (x / p["delta"]) ** 2)


def test_the_halo_relation_is_within_the_stellar_mass_systematics_of_two_independent_fits_from_1e12():
    """Stellar masses differ by 0.25 dex between surveys before the IMF (Leauthaud et al. 2012); the generator sits 0.04 to 0.26 dex
    under both from 1e12.5 up, 0.13 under Leauthaud and 0.12 over Behroozi at 1e12, so 0.30 dex is the line a later change must not cross."""
    virial_overdensity = 18 * math.pi ** 2 + 82 * (OMEGA_MATTER - 1.0) - 39 * (OMEGA_MATTER - 1.0) ** 2
    for log_halo in (12.0, 12.5, 13.0, 13.5, 14.0, 14.5):
        halo = 10 ** log_halo
        generator = math.log10(mass_chain.median_stellar_mass(halo))
        leauthaud = _leauthaud_log_stellar(math.log10(_mass_at_overdensity(halo, 200.0 * OMEGA_MATTER)))
        behroozi = _behroozi_log_stellar(math.log10(_mass_at_overdensity(halo, virial_overdensity)))
        assert abs(generator - leauthaud) <= 0.30, (log_halo, generator, leauthaud)
        assert abs(generator - behroozi) <= 0.30, (log_halo, generator, behroozi)
