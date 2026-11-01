# planetgen/physics/taxon_kernel.py

"""
The taxa of a belt from the disk's ice lines and the early heating of its bodies, in place of the Solar System's zone table (GEN.473, part 1;
docs/design/composition/small-bodies-composition.md section 13).

`small_body_composition.ZONES` is the one belt of the Solar System, in units of its snow line ([R1065: Table 6]): every taxon
row is a flagged stand-in. This module derives the same kind of kernel, as a function of the disk temperature a belt formed at and the
aluminium-26 its star was born with, and returns the measured zone mix as one realisation of it (a test, never an input):

1. **Where the families change** are the ice lines of [R2726] (`condensate_shares`): water ice condenses below 182 K, so bodies formed
   colder than that accrete water and alter to the hydrated carbonaceous taxa (C, B) instead of staying dry (S, V, A, M, K, L); ammonia hydrate
   below 131 K, where the primitive taxon P gives way to the ice-and-organic-rich D. The temperature at a distance is the blackbody one
   (`condensate_shares.disk_temperature_k`), so a belt's position is its own star's: nothing is in units of the Solar System's 2.7 AU. The line between
   hydrated (C) and primitive (P) bodies has no condensation line of its own and is a fitted temperature (`PRIMITIVE_BELOW_K`).
2. **Mixing.** Giant planets scatter bodies between zones; the share of each zone's mass exchanged with its neighbours is `IMPLANT_FRACTION`
   (fitted to the C share of the inner belt and the S share of the outer one in Table 6; a belt without a giant has none, a wiring-step choice).
3. **Melting.** A body of the dry family melts and differentiates into a core (the taxon M), a mantle (A) and a crust (V) when its aluminium-26 heats it
   to its solidus: the heat per kilogram is the ratio times the aluminium in the body times the decay energy (`elements`), and a body accreted at
   time `t` after the isotope's formation has the ratio decayed by `exp(-t / tau)`. With accretion times spread uniformly over `ACCRETION_SPREAD_MYR`, the melted
   share is the time a body can wait and still melt over the spread: `melted_share(ratio)`, which is 0.55 at the Solar System's 5.2e-5 and zero below
   `melting_ratio()` (about 1.1e-5), so the star's birth environment (GEN.476, `birth_radioisotopes`) decides how igneous its belt is.

Limits. The families are four, and the taxa inside each are shared out by fixed fractions (`DRY_SPLIT`, `IGNEOUS_SPLIT`, `HYDRATED_SPLIT`: the E taxon, an igneous
reduced one, has no source of its own here and is zero); the hot refractory K and L taxa have no temperature of their own; the lines are Lodders's at 10^-4 bar and
the disk is an irradiated blackbody. The Solar System's zones are test points: the families agree with Table 6 to about 0.15 (rms of the share by zone). Nothing in generation
calls this yet: the wiring step replaces `small_body_composition.kernel` with it, handing it the belt's temperature and the star's aluminium-26 ratio.
"""

import math

from planetgen.physics import condensate_shares, elements
from planetgen.util.checks import finite_domain

WATER_LINE_K = condensate_shares.ICE_CONDENSATION_K["H2O"]
AMMONIA_LINE_K = condensate_shares.ICE_CONDENSATION_K["NH3"]

PRIMITIVE_BELOW_K = 166.0
"""float: The temperature below which a body that accreted water stays primitive (the taxon P, ice and organics kept) and above which it is altered to the hydrated
carbonaceous taxa; fitted so the C to P change of Table 6 falls at its zones (x about 1.2). [Stand-in TK2, tier T6, tweak L3.]"""

WIDTH_LN_T = 0.05
"""float: The width of each change in the logarithm of the disk temperature (a width of 5 percent in temperature, 10 percent in distance), the softness of a line given a
disk that is not a single-temperature blackbody. [Stand-in TK2, tier T6, tweak L3.]"""

IMPLANT_FRACTION = 0.08
"""float: The share of the mass of each family that giant-planet scattering has moved to its neighbour families. [Stand-in TK3, tier T6, tweak L2.]"""

ALUMINIUM_MASS_FRACTION = 0.0100
"""float: Aluminium in a dry chondritic body, mass fraction (CI 0.0086 to H chondrite 0.0113, recalled). [Stand-in TK4, tier T4, tweak L0.]"""

MELT_ENERGY_J_PER_KG = 1.2e6
"""float: The heat that brings a chondritic body from 150 K to its solidus and melts a fifth of it: 900 J/kg/K times 1250 K and 0.2 times 4e5 J/kg of fusion
(textbook values; order of magnitude). [Stand-in TK4, tier T5, tweak L0.]"""

ACCRETION_SPREAD_MYR = 3.0
"""float: The spread of the accretion times of the planetesimals of one zone, after the isotope formed: the irons of the Solar System accreted under 0.5 Myr and the
chondrites over 2 Myr (recalled). [Stand-in TK5, tier T4, tweak L3.]"""

ALUMINIUM_26_MEAN_LIFE_MYR = 1.035
"""float: The mean life of aluminium-26 in Myr (`birth_radioisotopes.MEAN_LIFE_YR`)."""

DRY_SPLIT = {"S": 0.85, "K": 0.10, "L": 0.05}
"""dict: The undifferentiated dry bodies by taxon (ordinary against CO and CV chondrite analogues). [Stand-in TK6, tier T9, tweak L4.]"""

IGNEOUS_SPLIT = {"M": 0.25, "A": 0.65, "V": 0.10}
"""dict: A melted body's mass by taxon: the core (the metal share of a chondrite, 0.2 to 0.3), the mantle and the crust (a tenth by partial melting). [Stand-in TK6, tier T6, tweak L0.]"""

HYDRATED_SPLIT = {"C": 0.85, "B": 0.15}
"""dict: The hydrated carbonaceous taxa by letter (the C to B proportion of the middle and outer belt). [Stand-in TK6, tier T9, tweak L4.]"""

FAMILIES = ("dry", "hydrated", "primitive", "ice-rich")
TAXA = ("E", "C", "B", "S", "A", "L", "K", "V", "P", "M", "D")


@finite_domain()
def melting_ratio():
    """The aluminium-26 over aluminium-27 ratio at which a body accreted at once heats to its solidus: `Q / (X_Al / M_27 N_A E)`, about 1.1e-5."""
    isotope = elements.isotope("Al-26")
    atoms_per_kg = ALUMINIUM_MASS_FRACTION / elements.element("Al").molar_mass * elements.AVOGADRO     # molar mass is in kg/mol
    return MELT_ENERGY_J_PER_KG / (atoms_per_kg * isotope.decay_energy_j)


@finite_domain()
def melted_share(ratio):
    """The share of the dry-zone bodies that melted at an aluminium-26 ratio at the isotope's formation: the latest accretion time that still melts over the spread."""
    threshold = melting_ratio()
    if ratio <= threshold:
        return 0.0
    latest = ALUMINIUM_26_MEAN_LIFE_MYR * math.log(ratio / threshold)
    return min(latest / ACCRETION_SPREAD_MYR, 1.0)


def _above(temperature_k, line_k):
    """The weight of being warmer than a line: logistic in the logarithm of the temperature, 0.5 at the line."""
    z = (math.log(temperature_k) - math.log(line_k)) / WIDTH_LN_T
    return 1.0 / (1.0 + math.exp(-z)) if z > -700 else 0.0


@finite_domain()
def family_shares(temperature_k, implant=IMPLANT_FRACTION):
    """The mass shares of the four families (dry, hydrated, primitive, ice-rich) at a disk temperature, after the mixing between neighbours."""
    if not temperature_k > 0.0:
        raise ValueError("the temperature must be positive")
    dry = _above(temperature_k, WATER_LINE_K)
    wet = 1.0 - dry
    hydrated = wet * _above(temperature_k, PRIMITIVE_BELOW_K)
    cold = wet - hydrated
    primitive = cold * _above(temperature_k, AMMONIA_LINE_K)
    shares = {"dry": dry, "hydrated": hydrated, "primitive": primitive, "ice-rich": cold - primitive}
    mixed = {name: (1.0 - implant) * share for name, share in shares.items()}
    for i, name in enumerate(FAMILIES):
        neighbours = [FAMILIES[j] for j in (i - 1, i + 1) if 0 <= j < len(FAMILIES)]
        for other in neighbours:
            mixed[other] += implant * shares[name] / len(neighbours)
    return mixed


@finite_domain()
def kernel(temperature_k, aluminium_26_ratio, implant=IMPLANT_FRACTION):
    """
    The taxon weights of a belt formed at this disk temperature round a star born with this aluminium-26 over aluminium-27 ratio, adding to 1.

    The dry family splits into undifferentiated taxa and the melted share's core, mantle and crust; the hydrated into C and B; P and D are one taxon each.
    """
    families = family_shares(temperature_k, implant)
    melted = melted_share(aluminium_26_ratio)
    weights = dict.fromkeys(TAXA, 0.0)
    for taxon, share in DRY_SPLIT.items():
        weights[taxon] += families["dry"] * (1.0 - melted) * share
    for taxon, share in IGNEOUS_SPLIT.items():
        weights[taxon] += families["dry"] * melted * share
    for taxon, share in HYDRATED_SPLIT.items():
        weights[taxon] += families["hydrated"] * share
    weights["P"] += families["primitive"]
    weights["D"] += families["ice-rich"]
    return weights


@finite_domain()
def kernel_at_distance(distance_au, luminosity_w, aluminium_26_ratio, implant=IMPLANT_FRACTION):
    """The kernel of a belt at a distance from a star of this luminosity (the blackbody disk temperature there)."""
    return kernel(condensate_shares.disk_temperature_k(distance_au, luminosity_w), aluminium_26_ratio, implant)
