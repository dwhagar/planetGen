# planetgen/physics/planetary_nebula_ejecta.py

"""
What a planetary nebula is made of (GEN.1270, part b of GEN.1268, reports/gen1268-nebula-composition-2026-11-05.md section 5): the progenitor star
that lost it, and the helium, carbon, nitrogen and oxygen of its ejecta after the dredge-ups.

A planetary nebula stores no progenitor, so one is drawn. The draw is **one stream of the nebula's own key** (`progenitor_mass_sol`): the same
nebula always has the same progenitor, and a central star or white dwarf later read from the same key and class gets the same mass. The values
are for the composition item (GEN.1276) to write, with the progenitor mass, so they are returned together (`PlanetaryNebulaEjecta`).

**The progenitor.** An initial mass drawn from the initial mass function (`imf.sample`) between `MIN_MASS_SOL` and the top of the white dwarf
relation, inside the window the class allows:

* H (young) and J (evolved): any mass.
* K (carbon-rich): from the least mass that dredges carbon up at the place's metallicity (`agb_carbon.tdu_minimum_mass_sol`) to where hot
  bottom burning (`agb_carbon.HBB_MIN_MASS_SOL`) takes the carbon again, and the ejecta must come out with C/O at or above 1 (the draw is
  repeated, from the same stream, until they do; a fixed number of tries, then the heaviest mass of the window).
* L (nitrogen-rich bipolar, Type I): from `HBB_MIN_MASS_SOL` up, where hot bottom burning makes nitrogen.

**The birth composition** is the thin-disk mean [X/Fe] at the formation-epoch [Fe/H]: the place's [Fe/H] now, less the drift of the
disk's gas over the progenitor's life (`age_metallicity.THIN_SLOPE_DEX_PER_GY` times the lifespan from the star track), so an older progenitor
formed from poorer gas.

**The ejecta**, in this order, as mass fractions (every step keeps the sum at 1):

1. **First dredge-up** (`stellar_abundances.after_first_dredge_up`): helium up, carbon down, nitrogen up. The table ends at 2.5 solar masses; a
   heavier star takes its edge, as the star sheet does.
2. **Second dredge-up** (a star of `SECOND_DREDGE_UP_TABLE`'s lightest row or more): helium up by `second_dredge_up_delta_y` (the table of
   [R5098: Table 1], linear between its rows), taken from hydrogen.
3. **Third dredge-up** (a star above `agb_carbon.tdu_minimum_mass_sol`): the carbon, helium and oxygen of the intershell
   (`INTERSHELL_*`) brought up over the thermal pulses (`agb_carbon.dredged_core_mass_sol`), all of which ends in the ejecta, so the mass dredged
   `M_dred` is mixed with the envelope the thermally pulsing branch loses (`M_ej`, the mass after the red giant branch less the white dwarf):
   each fraction is `(1 - M_dred/M_ej)` of the envelope's plus `M_dred/M_ej` of the intershell's. This is the mass-weighted average of everything lost on the branch;
   the nebula itself is the late part, richer in the dredged carbon (the first build does not separate them).
4. **Hot bottom burning** (a star of `HBB_MIN_MASS_SOL` or more): `HBB_CARBON_TO_NITROGEN` of the carbon becomes nitrogen (14.007/12.011 as much
   mass, the difference taken from hydrogen).

Neon, sulfur, argon, iron and the rest stay at the birth composition, which is how planetary nebula abundances are used observationally.
**Type I** is the usual definition, helium over hydrogen above `TYPE_I_HE_TO_H` or log N/O above `TYPE_I_LOG_N_TO_O`; the result flags it.

**The white dwarf** is the initial-final relation of `remnant_mass`; the ejected mass is the initial mass less it.
"""

import math
from dataclasses import dataclass

from planetgen import tuning
from planetgen.physics import age_metallicity, agb_carbon, imf, post_main_sequence, remnant_mass
from planetgen.physics import stellar_abundances as sa
from planetgen.util import draw

H, J, K, L = "H", "J", "K", "L"
CLASSES = (H, J, K, L)
"""tuple: The nebula classes this module covers: H young and J evolved planetary nebulae, K carbon-rich, L nitrogen-rich bipolar."""

MIN_MASS_SOL = 0.9
"""float: The lightest progenitor drawn. [Stand-in PNE1, tier T6, tweak L1: a star must end its life within the age of the disk, about 10 Gy,
which is a little under a solar mass at the track's lifespans; the initial-final relation starts at 0.83.]"""

MAX_MASS_SOL = tuning.WD_IFMR_PIECES[-1][1]
"""float: The heaviest progenitor drawn: the top of the initial-final relation of white dwarfs (7.2 solar masses), where core collapse begins."""

RGB_MASS_LOSS_SOL = 0.2
"""float: Mass lost on the red giant branch before the thermally pulsing branch begins. [Stand-in PNE2, tier T4, tweak L0: [R5098] says only that the
loss on the branch may be "as much as 30 percent of the star's total mass for the lowest mass stars", which brackets 0.2 for a star of 1 to 2 solar masses
(0.3 to 0.6) from above; the per-mass numbers are in the on-line tables of [R5100], which could not be read, so the value stays the recalled one.]"""

SECOND_DREDGE_UP_TABLE = ((4.0, 0.0), (4.5, 0.004), (5.0, 0.018), (5.5, 0.033), (6.0, 0.044), (7.0, 0.059), (8.0, 0.066))
"""tuple: `(initial mass, solar masses; helium mass fraction the second dredge-up adds)`, the rise in surface Y from after the first to after the second
dredge-up in [R5098: Table 1] (Z = 0.02: Y 0.293 to 0.297 at 4.5, 0.291 to 0.309 at 5.0, 0.289 to 0.322 at 5.5, 0.289 to
0.333 at 6.0, 0.291 to 0.350 at 7.0, 0.296 to 0.362 at 8.0; no change from 1 to 4.0). [Verified, tier T1: the table is for solar metallicity
only and the same rise is used at every metallicity, so the lower-Z shift is not modelled; above 7.2 solar masses nothing is drawn.]"""

INTERSHELL_HE = 0.75
INTERSHELL_O = 0.02
INTERSHELL_REST = 0.01
"""float: The helium, oxygen and remaining (taken as neon-22) mass fractions of the intershell dredged up; the carbon is
`agb_carbon.INTERSHELL_CARBON_FRACTION`, and the four add to 1. [Stand-in PNE4, tier T4, tweak L0: recalled, intershell 70 to 75 percent helium, 20
to 25 percent carbon, about 2 percent oxygen.]"""

HBB_CARBON_TO_NITROGEN = 0.9
"""float: The share of the envelope's carbon hot bottom burning turns into nitrogen. [Stand-in PNE5, tier T6, tweak L3: set so that a solar-metallicity star of 5 solar
masses comes out Type I (log N/O above -0.3), as the nitrogen-rich bipolar class is defined.]"""

TYPE_I_HE_TO_H = 0.125
TYPE_I_LOG_N_TO_O = -0.3
"""float: The limits of Type I planetary nebulae: He/H (by number) above the first or log N/O above the second. [Verified, tier T1: the
classification as quoted in [R5099]: He/H of 0.125 or more, or N/O above 0.5 (log N/O above -0.30); the
Magellanic Clouds and another author (N/O above 0.8) use other limits, which are not used here.]"""

CARBON_RICH_TRIES = 40
"""int: Draws made from the stream for a class K progenitor before the heaviest mass of the window is taken."""

_NITROGEN_PER_CARBON = 14.007 / 12.011


@dataclass(frozen=True)
class PlanetaryNebulaEjecta:
    """
    A planetary nebula's progenitor and ejecta.

    Attributes:
        nebula_class (str): `H`, `J`, `K` or `L`.
        progenitor_mass_sol (float): The progenitor's initial mass.
        white_dwarf_mass_sol (float): The white dwarf it leaves.
        ejected_mass_sol (float): The mass it lost: the initial mass less the white dwarf.
        birth_fe_h (float): [Fe/H] of the gas it formed from.
        pulses (int): Thermal pulses over its asymptotic branch (0 when it does not dredge carbon up).
        dredged_mass_sol (float): The intershell mass dredged up over them.
        third_dredge_up, hot_bottom_burning (bool): Whether each ran.
        composition (StellarComposition): The ejecta's mass fractions (`x`, `y`, `z`, and `elements`).
        c_to_o, n_to_o (float): Number ratios.
        he_to_h (float): Helium over hydrogen by number.
        carbon_rich, type_i (bool): C/O at or above 1; the Type I definition.
    """

    nebula_class: str
    progenitor_mass_sol: float
    white_dwarf_mass_sol: float
    ejected_mass_sol: float
    birth_fe_h: float
    pulses: int
    dredged_mass_sol: float
    third_dredge_up: bool
    hot_bottom_burning: bool
    composition: object
    c_to_o: float
    n_to_o: float
    he_to_h: float
    carbon_rich: bool
    type_i: bool

    def values(self):
        """The values the composition item stores, as a dict: the progenitor and white dwarf masses, the helium, hydrogen and metal mass fractions, the
        mass fractions of C, N, O, Ne, S, Ar and Fe, and the ratios."""
        c = self.composition
        return {"progenitor_mass_sol": self.progenitor_mass_sol, "white_dwarf_mass_sol": self.white_dwarf_mass_sol, "ejected_mass_sol": self.ejected_mass_sol,
                "mass_fraction_h": c.x, "mass_fraction_he": c.y, "mass_fraction_z": c.z,
                **{f"mass_fraction_{e.lower()}": c.elements[e] for e in ("C", "N", "O", "Ne", "S", "Ar", "Fe")},
                "c_to_o": self.c_to_o, "n_to_o": self.n_to_o, "he_to_h": self.he_to_h, "carbon_rich": self.carbon_rich, "type_i": self.type_i}


def _stream(key):
    return draw.Stream(f"planetary-nebula-progenitor|{key}")


def lifespan_gy(initial_mass_sol, fe_h):
    """The age, Gy, at which a star of this mass and [Fe/H] ends its life (the end of its track)."""
    return post_main_sequence.timeline(initial_mass_sol, fe_h).end_gy


def birth_fe_h(place_fe_h, initial_mass_sol):
    """[Fe/H] of the gas a progenitor formed from: the place's now, less the disk's drift over the progenitor's lifespan."""
    return place_fe_h + age_metallicity.THIN_SLOPE_DEX_PER_GY * lifespan_gy(initial_mass_sol, place_fe_h)


def second_dredge_up_delta_y(initial_mass_sol):
    """The helium mass fraction the second dredge-up adds at an initial mass: `SECOND_DREDGE_UP_TABLE`, linear between its rows, 0 at and under its
    lightest row, and its heaviest row's value above."""
    table = SECOND_DREDGE_UP_TABLE
    if initial_mass_sol <= table[0][0]:
        return 0.0
    for (m0, y0), (m1, y1) in zip(table, table[1:]):
        if initial_mass_sol <= m1:
            return y0 + (y1 - y0) * (initial_mass_sol - m0) / (m1 - m0)
    return table[-1][1]


def mass_window_sol(nebula_class, place_fe_h):
    """The `(least, greatest)` initial mass a class allows (see the module doc)."""
    if nebula_class not in CLASSES:
        raise ValueError(f"unknown planetary nebula class {nebula_class!r}; expected one of {', '.join(CLASSES)}")
    if nebula_class == K:
        return max(MIN_MASS_SOL, agb_carbon.tdu_minimum_mass_sol(place_fe_h)), agb_carbon.HBB_MIN_MASS_SOL
    if nebula_class == L:
        return agb_carbon.HBB_MIN_MASS_SOL, MAX_MASS_SOL
    return MIN_MASS_SOL, MAX_MASS_SOL


def ejecta_of(nebula_class, initial_mass_sol, place_fe_h):
    """The `PlanetaryNebulaEjecta` of a progenitor of a given mass (no draw): the dredge-ups of the module doc."""
    fe_h = birth_fe_h(place_fe_h, initial_mass_sol)
    star = sa.composition(fe_h, sa.abundance_offsets(fe_h, sa.BRANCH_THIN, (0.0,) * len(sa.RESIDUAL_ELEMENTS)))
    envelope = sa.after_first_dredge_up(star, initial_mass_sol)
    elements = dict(envelope.elements)
    white_dwarf = remnant_mass.white_dwarf_mass_sol(initial_mass_sol)
    hbb = initial_mass_sol >= agb_carbon.HBB_MIN_MASS_SOL
    delta_y = second_dredge_up_delta_y(initial_mass_sol)
    elements["He"] += delta_y
    elements["H"] -= delta_y
    pulses, dredged = 0, 0.0
    third = initial_mass_sol >= agb_carbon.tdu_minimum_mass_sol(fe_h)
    if third:
        pulses = int(agb_carbon.pulses_in_branch(initial_mass_sol))
        dredged = agb_carbon.dredged_core_mass_sol(initial_mass_sol, pulses)
        lost = max(initial_mass_sol - RGB_MASS_LOSS_SOL - white_dwarf, dredged)
        share = dredged / lost
        carbon = agb_carbon.INTERSHELL_CARBON_FRACTION
        rest = 1.0 - carbon - INTERSHELL_HE - INTERSHELL_O
        assert abs(rest - INTERSHELL_REST) < 1e-9, "the intershell fractions must add up to 1"
        for symbol in elements:
            elements[symbol] *= 1.0 - share
        elements["He"] += share * INTERSHELL_HE
        elements["C"] += share * carbon
        elements["O"] += share * INTERSHELL_O
        elements["Ne"] += share * rest
    if hbb:
        burned = HBB_CARBON_TO_NITROGEN * elements["C"]
        elements["C"] -= burned
        elements["N"] += burned * _NITROGEN_PER_CARBON
        elements["H"] -= burned * (_NITROGEN_PER_CARBON - 1.0)
    metals = sum(value for symbol, value in elements.items() if symbol not in ("H", "He"))
    composition = sa.StellarComposition(fe_h, elements["H"], elements["He"], metals, elements, envelope.offsets)
    c_to_o = (elements["C"] / 12.011) / (elements["O"] / 15.999)
    n_to_o = (elements["N"] / 14.007) / (elements["O"] / 15.999)
    he_to_h = (elements["He"] / 4.0026) / (elements["H"] / 1.008)
    return PlanetaryNebulaEjecta(
        nebula_class, initial_mass_sol, white_dwarf, initial_mass_sol - white_dwarf, fe_h, pulses, dredged, third, hbb, composition, c_to_o, n_to_o,
        he_to_h, c_to_o >= 1.0, he_to_h > TYPE_I_HE_TO_H or (n_to_o > 0.0 and math.log10(n_to_o) > TYPE_I_LOG_N_TO_O))


def progenitor(nebula_class, place_fe_h, key):
    """
    The `PlanetaryNebulaEjecta` of the nebula `key` (its object ID) of class `nebula_class` at a place of [Fe/H] `place_fe_h`: the progenitor drawn
    from the stream of the key, inside the class's window, and its ejecta. The same arguments always give the same answer; a central star or a white
    dwarf that reads the same key, class and place gets the same progenitor mass.
    """
    low, high = mass_window_sol(nebula_class, place_fe_h)
    stream = _stream(key)
    if nebula_class != K:
        return ejecta_of(nebula_class, imf.sample(low, high, rng=stream), place_fe_h)
    for _ in range(CARBON_RICH_TRIES):
        found = ejecta_of(nebula_class, imf.sample(low, high, rng=stream), place_fe_h)
        if found.carbon_rich:
            return found
    return ejecta_of(nebula_class, 0.99 * high, place_fe_h)


def progenitor_mass_sol(nebula_class, place_fe_h, key):
    """The initial mass of the progenitor of the nebula `key` (see `progenitor`)."""
    return progenitor(nebula_class, place_fe_h, key).progenitor_mass_sol
