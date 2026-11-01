# planetgen/physics/hii.py

"""
H II regions from their stars' ionizing photons (GEN.150).

The rules are docs/design/systems/nebula-and-asteroid-field-classes.md, "H II radius
from the ionizing photon rate" and "Lifetimes by class":

1. **Photons.** Each hot main-sequence star gives off `Q` photons a second
   above 13.6 eV (`log_q_of_type`, `log_q_of_mass`): the recalled
   [R2915] values for O3V to B5V, which fall by
   about a decade per subclass from B0 to B2. Stars under 5.5 Msun and
   cooler than B5 ionize nothing (`Q` 0).
2. **Radius.** The Stromgren sphere, `R_S = 3.1 pc (Q / 1e49)^(1/3) (n /
   100)^(-2/3)`, then the front keeps growing as a D-type front (Spitzer),
   `R(t) = R_S (1 + 7 c_s t / (4 R_S))^(4/7)` with `c_s` 10 km/s, while the
   density behind it falls as `n (R_S / R)^(3/2)`. `region` draws nothing: it
   returns the radius and density of a cloud of density `n` ionized by `Q`
   for `t`.
3. **Clusters.** A cloud's `Q` comes from its stars: `cluster_hosts` draws a
   cluster of a given mass from the Kroupa IMF (`tuning.IMF_*`), keeps the
   stars still on the main sequence at the cloud's age and adds their `Q`.
   About 5.7e46 photons a second per solar mass of cluster, against the
   design note's 5.4e46.
"""

import math

from planetgen.physics import imf, stellar_evolution
from planetgen.util import draw
from planetgen.physics import unit_registry as units

_pc_to_ly = units.converter("pc", "ly")

ALPHA_B_CM3_S = 2.6e-13
"""float: Case B recombination coefficient of hydrogen at 1e4 K."""

SOUND_SPEED_PC_MYR = 10.0 * 1.0227
"""float: The ionized gas's sound speed, 10 km/s, in pc per Myr (1 km/s is 1.0227 pc/Myr)."""

PC_CM = 3.0857e18
LY_PER_PC = _pc_to_ly(1.0)
"""float: Light-years in one parsec, exact from the two definitions in `physics/constants.py` (GEN.572)."""

MIN_IONIZING_MASS_SOL = 5.5
"""float: The lightest star counted as a host: a B5V, Q about 3e43 a second. Cooler
stars light dust, they do not ionize it (Q is 0 and they are no host)."""

TYPE_ANCHORS = (
    ("O3V", 60.0, 49.64), ("O5V", 40.0, 49.22), ("O7V", 25.0, 48.8), ("O9V", 18.0, 48.2),
    ("B0V", 15.0, 47.6), ("B1V", 11.0, 46.8), ("B2V", 8.5, 45.7), ("B3V", 7.5, 45.0), ("B5V", 5.5, 43.5),
)
"""tuple: `(type, zero-age mass in Msun, log10 Q)` from the hottest down: O3V to B5V
(design note, "H II radius from the ionizing photon rate"; the masses and the
log Q values are recalled to about 0.15 dex for O stars and 0.3 for B stars).
Above O3V the photon rate stops rising (a plateau at the O3V value)."""

_SUBCLASS_STEPS = {"O": 0, "B": 10}
"""dict: The position of a letter's subclass 0 on the O3V-to-B5V ladder (O3 is 3, B0 is 10)."""


def _ladder_index(letter, subclass):
    """A spectral type as a number: O3 is 3, O9 is 9, B0 is 10, B5 is 15."""
    return _SUBCLASS_STEPS[letter] + subclass


_ANCHOR_INDEX = (3, 5, 7, 9, 10, 11, 12, 13, 15)


def log_q_of_type(letter, subclass):
    """
    log10 of the ionizing photons a second of a main-sequence star of
    spectral type `letter` + `subclass` ("O", 5 is O5V), linear in the
    subclass between the anchors; `None` for a type that ionizes nothing
    (B6 and later, A and cooler).
    """
    if letter not in _SUBCLASS_STEPS:
        return None
    index = _ladder_index(letter, subclass)
    if index > _ANCHOR_INDEX[-1]:
        return None
    if index <= _ANCHOR_INDEX[0]:
        return TYPE_ANCHORS[0][2]
    for (i0, (_, _, q0)), (i1, (_, _, q1)) in zip(zip(_ANCHOR_INDEX, TYPE_ANCHORS),
                                                  zip(_ANCHOR_INDEX[1:], TYPE_ANCHORS[1:])):
        if index <= i1:
            return q0 + (q1 - q0) * (index - i0) / (i1 - i0)
    return None


def log_q_of_mass(mass_sol):
    """log10 Q for a main-sequence star of zero-age mass `mass_sol`: log-linear between
    the anchors, flat above O3V's 60 Msun, `None` under `MIN_IONIZING_MASS_SOL`."""
    if mass_sol < MIN_IONIZING_MASS_SOL:
        return None
    if mass_sol >= TYPE_ANCHORS[0][1]:
        return TYPE_ANCHORS[0][2]
    for (_, m_high, q_high), (_, m_low, q_low) in zip(TYPE_ANCHORS, TYPE_ANCHORS[1:]):
        if mass_sol >= m_low:
            fraction = (math.log(mass_sol) - math.log(m_low)) / (math.log(m_high) - math.log(m_low))
            return q_low + fraction * (q_high - q_low)
    return TYPE_ANCHORS[-1][2]


def type_of_mass(mass_sol):
    """The anchor type (`"O5V"`) nearest a star's mass, in log mass; `None` for a star
    under `MIN_IONIZING_MASS_SOL`."""
    if mass_sol < MIN_IONIZING_MASS_SOL:
        return None
    return min(TYPE_ANCHORS, key=lambda anchor: abs(math.log(mass_sol) - math.log(anchor[1])))[0]


def mass_range_of_type(name, upper_sol):
    """
    `(low, high)` solar masses of the stars `type_of_mass` names `name`: from the log-midpoint to the next lighter anchor
    to the one to the next heavier. The hottest anchor reaches `upper_sol` and the coolest `MIN_IONIZING_MASS_SOL`.

    Raises:
        ValueError: `name` is not an anchor type.
    """
    names = [anchor[0] for anchor in TYPE_ANCHORS]
    if name not in names:
        raise ValueError(f"{name!r} is not one of {names}")
    index = names.index(name)
    mass = TYPE_ANCHORS[index][1]
    high = upper_sol if index == 0 else math.sqrt(mass * TYPE_ANCHORS[index - 1][1])
    low = MIN_IONIZING_MASS_SOL if index == len(TYPE_ANCHORS) - 1 else math.sqrt(mass * TYPE_ANCHORS[index + 1][1])
    return low, high


def type_name_of_mass(mass_sol):
    """A main-sequence star's spectral type from its mass, to the nearest anchor type
    (`type_of_mass`) down to B5V, then B8V under 5.5 Msun and A0V under 3.8 Msun."""
    named = type_of_mass(mass_sol)
    if named is not None:
        return named
    return "B8V" if mass_sol >= 3.8 else "A0V"


def stromgren_radius_pc(q, density_cm3):
    """`R_S`, parsecs, of a uniform cloud of hydrogen density `density_cm3` ionized by `q`
    photons a second: `(3 Q / (4 pi n^2 alpha_B))^(1/3)`."""
    if q <= 0.0 or density_cm3 <= 0.0:
        return 0.0
    return (3.0 * q / (4.0 * math.pi * density_cm3 ** 2 * ALPHA_B_CM3_S)) ** (1.0 / 3.0) / PC_CM


def expanded_radius_pc(stromgren_pc, age_myr):
    """The radius after the front has expanded for `age_myr` (Spitzer): `R_S (1 + 7 c_s t / (4 R_S))^(4/7)`."""
    if stromgren_pc <= 0.0:
        return 0.0
    return stromgren_pc * (1.0 + 7.0 * SOUND_SPEED_PC_MYR * max(age_myr, 0.0) / (4.0 * stromgren_pc)) ** (4.0 / 7.0)


def region(q, density_cm3, age_myr):
    """
    The H II region ionized by `q` photons a second in a cloud of density
    `density_cm3` for `age_myr`: `(radius_ly, density_cm3)`, the radius
    `expanded_radius_pc` and the density `n (R_S / R)^(3/2)` behind the front.
    `(0.0, density_cm3)` when nothing ionizes.
    """
    r_s = stromgren_radius_pc(q, density_cm3)
    if r_s <= 0.0:
        return 0.0, density_cm3
    radius = expanded_radius_pc(r_s, age_myr)
    return radius * LY_PER_PC, density_cm3 * (r_s / radius) ** 1.5


def _poisson(mean):
    """A Poisson draw on the run's stream: Knuth's product for a small mean, the normal
    approximation (continuity-corrected) above 30."""
    if mean <= 0.0:
        return 0
    if mean > 30.0:
        return max(0, int(round(draw.gauss(mean, math.sqrt(mean)))))
    limit, product, count = math.exp(-mean), draw.random(), 0
    while product > limit:
        count += 1
        product *= draw.random()
    return count


def _stars_per_solar_mass():
    """The mean number of host-capable stars (at least `MIN_IONIZING_MASS_SOL`) per
    solar mass of cluster: the IMF's count above that mass over its total mass."""
    low, high = imf.lowest_mass_sol(), imf.highest_mass_sol()
    return imf.integral(0.0, MIN_IONIZING_MASS_SOL, high) / imf.integral(1.0, low, high)


def cluster_hosts(cluster_mass_sol, age_myr):
    """
    The ionizing stars of a cluster of `cluster_mass_sol` solar masses at
    `age_myr`: the number of stars of at least `MIN_IONIZING_MASS_SOL` is
    Poisson with the IMF's mean for that cluster mass, each draws a mass from
    the IMF, and the ones whose main-sequence life is over by `age_myr`
    (they are supernovae by now) are dropped.

    Returns:
        list: The surviving stars' masses in solar masses, heaviest first.
    """
    count = _poisson(cluster_mass_sol * _stars_per_solar_mass())
    age_gy = age_myr / 1000.0
    masses = [imf.sample(MIN_IONIZING_MASS_SOL) for _ in range(count)]
    return sorted((m for m in masses if stellar_evolution.main_sequence_lifetime_gy(m) > age_gy), reverse=True)


def total_log_q(masses, extra_log_q=()):
    """log10 of the summed photon rate of stars (`masses`, solar masses) and of any
    `extra_log_q` values; `None` when nothing ionizes."""
    rates = [10.0 ** log_q_of_mass(m) for m in masses if log_q_of_mass(m) is not None]
    rates.extend(10.0 ** value for value in extra_log_q)
    return math.log10(sum(rates)) if rates else None


def host_summary(counts):
    """`{type: count}` as the stored text, hottest first: `"O5V x1, O9V x4, B0V x9"`."""
    order = [anchor[0] for anchor in TYPE_ANCHORS]
    ranked = sorted(counts.items(), key=lambda item: order.index(item[0]) if item[0] in order else len(order))
    return ", ".join(f"{name} x{count}" for name, count in ranked if count)


def host_counts(masses):
    """`{anchor type: count}` for stars of `masses`."""
    counts = {}
    for mass in masses:
        name = type_of_mass(mass)
        if name is not None:
            counts[name] = counts.get(name, 0) + 1
    return counts
