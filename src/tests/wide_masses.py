# tests/wide_masses.py

"""Stellar masses over the wide range (3e9 to 2e11 solar masses, the galaxy mass function's weight) for tests that
check how a law changes with mass: the Auto window (`mass_chain.AUTO_LOG_MASS_RANGE`) is a factor of two."""

import math

from planetgen.galaxy import mass_chain

WIDE_LOG_MASS_RANGE = (math.log10(3.0e9), math.log10(2.0e11))


def wide_mass(seed):
    """The stellar mass, solar masses, `seed` would have drawn before the Auto window narrowed (GEN.603)."""
    share = mass_chain._stream(seed, "stellar-mass").random()
    return 10.0 ** mass_chain.log_stellar_mass_at(share, WIDE_LOG_MASS_RANGE)
