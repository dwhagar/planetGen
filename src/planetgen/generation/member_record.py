# planetgen/generation/member_record.py

"""
The member record: one evolution of a star, frozen (GEN.1009, GEN.1006 section 4.1).

A star-like member of a system used to be evolved twice, once by the tree that draws the orbits (`multistar._star_leaf`) and again
by the build that makes its `Star`, each with its own random source and the build with the primary's [Fe/H] where the tree had the
reference one. The two could disagree on the radius by several per cent, enough to put a pair drawn at its Roche floor inside a
star. `make` is the one place a birth mass, an age and a [Fe/H] become a state, under a stream of the member's own seed, so the
state is a function of those four and nothing else: the tree reads the record, the build makes its star from the same record
(`Star.from_record`), and a stored tree, which keeps the seed and the birth mass, rebuilds the same record at any later time.

A record carries what the rest of the generator needs from the state: the present mass and radius, the largest radius and mass the
member has had (`radius_max_au`, `mass_max_solar`: the red-giant and AGB tips of a star that is smaller now, for the history rule),
and for a star the full `star_params` the build makes its star from.
"""

import dataclasses

from planetgen.physics import (alpha_tracks, brown_dwarf, constants, helium_core_white_dwarf, post_main_sequence, pre_main_sequence,
                              stellar_abundances)
from planetgen.physics import stellar_evolution as evolution
from planetgen.physics import unit_registry as units
from planetgen.util import draw

_km_to_m = units.converter("km", "m")
_m_to_au = units.converter("m", "au")
_kg_to_msun = units.converter("kg", "Msun")

LIGHTEN_FACTOR = 0.8
"""float: A birth mass that has already collapsed at the system's age is cut by this factor until it is alive (the same cut
`Star` makes for a companion of a specified-type primary, GEN.53)."""

YOUNG_KINDS = ("main-sequence", "pre-main-sequence", "protostar")
"""tuple: The kinds that have not left the main sequence: their largest radius is the present one."""

STAR_KINDS = ("main-sequence", "pre-main-sequence", "protostar", "giant", "white-dwarf")
"""tuple: The kinds `make` gives a star (a `Star` can be made from the record); a brown dwarf is the other thing it makes."""


@dataclasses.dataclass(frozen=True)
class MemberRecord:
    """
    One member's evolved state at the system's age. `kind` is `brown-dwarf`, `protostar`, `pre-main-sequence`, `main-sequence`,
    `giant`, `white-dwarf`, or (from `compact`) `neutron-star` or `black-hole`; `state` the Yerkes class of a star
    (`V`, `III`, `VII`, ...), `PMS` or `protostar` for a young one, the spectral type of a brown dwarf, and None for a remnant.
    `mass_solar` and `radius_au` are the present values and `initial_mass_solar` the birth mass the evolution took (a lightened one
    if the first had collapsed). `radius_max_au` is the largest radius the member has had after its main sequence (the present one
    for a member that has not left it), `mass_max_solar` the largest mass it has had. `params` is the `star_params` dict a star is
    made from, as a tuple of items so the record stays hashable; None for a brown dwarf or a compact remnant.
    """
    kind: str
    state: str | None
    mass_solar: float
    initial_mass_solar: float
    radius_au: float
    radius_max_au: float
    mass_max_solar: float
    teff_k: float | None
    luminosity_sol: float | None
    age_yr: float
    fe_h: float
    seed: int
    params: tuple | None = None

    def star_params(self):
        """The `star_params` dict a `Star` is made from, `None` for a record that is not a star."""
        return None if self.params is None else dict(self.params)


def stream(seed):
    """The random stream of the member with this `seed`: everything random in its evolution (a horizontal-branch colour, a
    Wolf-Rayet luminosity) is drawn from it, so the state is the same wherever it is asked for."""
    return draw.Stream(f"member-record:{int(seed)}")


def _radius_au_of_km(radius_km):
    return _m_to_au(_km_to_m(radius_km))


def _kind_of(yerkes):
    """The member kind of a star of this Yerkes class (a contracting star's is its stage)."""
    if yerkes == pre_main_sequence.PROTOSTAR:
        return "protostar"
    if yerkes in pre_main_sequence.STAGE_CLASSES:
        return "pre-main-sequence"
    return "main-sequence" if yerkes in ("V", "VI") else "white-dwarf" if yerkes == "VII" else "giant"


def _brown_dwarf(initial_mass_solar, age_yr, fe_h, seed):
    """The record of a substellar member: the L, T or Y type its cooling gives at the age (GEN.252), and the radius the same laws fix."""
    dwarf = brown_dwarf.state(initial_mass_solar, age_yr / 1.0e9)
    radius = _m_to_au(dwarf["radius_sol"] * constants.SOLAR_RADIUS_M)
    return MemberRecord("brown-dwarf", dwarf["spectral_type"], initial_mass_solar, initial_mass_solar, radius, radius,
                        initial_mass_solar, dwarf["temperature_k"], dwarf["luminosity_sol"], age_yr, fe_h, seed)


REMNANT_AGE_YR = 1.0e13
"""float: An age at which every window of any star's track is finished, for the tips a remnant's progenitor passed. Source tier T6, tweak level L0."""


def _largest_radius_au(initial_mass_solar, age_yr, fe_h, present_radius_au):
    """The largest radius a star of this birth mass has had by `age_yr`: the larger of its present radius and its giant-branch tips."""
    try:
        line = post_main_sequence.timeline(initial_mass_solar, alpha_tracks.track_fe_h(fe_h, age_yr / 1.0e9))
    except ValueError:
        return present_radius_au
    tip_sol = line.radius_max_sol(age_yr / 1.0e9)
    return max(present_radius_au, _m_to_au(tip_sol * constants.SOLAR_RADIUS_M))


def progenitor_tip_au(initial_mass_solar, fe_h=0.0):
    """The largest radius, AU, a star born at this mass reaches before it becomes a remnant: its giant-branch tips, every window finished (0 for a mass that has no track)."""
    return _largest_radius_au(initial_mass_solar, REMNANT_AGE_YR, fe_h, 0.0)


def make(seed, initial_mass_solar, age_yr, fe_h, lighten=True, substellar=False, white_dwarf_mass_sol=None):
    """
    The record of the member born at `initial_mass_solar` and seen at `age_yr` with birth [Fe/H] `fe_h`, evolved under the stream of
    `seed`. A brown dwarf is below the hydrogen-burning limit. A star is `evolution.evolve_star`'s state (a contracting one's too);
    its radius is the one `star_params` gives, so the star made from the record has exactly the record's radius. A birth mass that has
    collapsed by the age is cut by `LIGHTEN_FACTOR` until it is alive when `lighten`, and raises `ValueError` when not. A `substellar`
    member was drawn from the brown dwarf mass function and stays one at any mass (the objects of 0.074 to 0.08 solar masses among them
    are the stable late-M and L dwarfs of the census); the limit is the one of the member's [Fe/H] (GEN.936). A star that is a white
    dwarf by the age is one of `white_dwarf_mass_sol` solar masses when given (a helium core, GEN.996); any other star ignores it.
    """
    if substellar or initial_mass_solar < stellar_abundances.hydrogen_burning_limit_for_fe_h(fe_h):
        return _brown_dwarf(initial_mass_solar, age_yr, fe_h, seed)
    rng = stream(seed)
    mass = initial_mass_solar
    while True:
        state = evolution.evolve_star(mass, age_yr / 1.0e9, rng, fe_h=fe_h, white_dwarf_mass_sol=white_dwarf_mass_sol)
        if state is not None or not lighten:
            break
        draw.redrawn("star.mass_lightened")
        mass *= LIGHTEN_FACTOR
    if state is None:
        raise ValueError(f"a star of {initial_mass_solar:g} solar masses has collapsed by {age_yr:g} years")
    params = evolution.star_params(mass, age_yr / 1.0e9, state)
    yerkes = params["yerkes_class"]
    kind = _kind_of(yerkes)
    radius = _radius_au_of_km(params["radius_km"])
    largest = radius if kind in YOUNG_KINDS else _largest_radius_au(mass, age_yr, fe_h, radius)
    if white_dwarf_mass_sol is not None and kind == "white-dwarf":
        # A helium core: the giant gave up its envelope at the lobe it filled, on the red giant branch, and never grew past it.
        largest = max(radius, _m_to_au(helium_core_white_dwarf.giant_radius_rsun(white_dwarf_mass_sol, fe_h)
                                       * constants.SOLAR_RADIUS_M))
    return MemberRecord(kind, yerkes, state["mass_sol"], mass, radius, largest, mass, params["temperature_k"],
                        state["luminosity_sol"], age_yr, fe_h, seed, tuple(sorted(params.items())))


def compact(kind, seed, mass_solar, progenitor_solar, radius_au, age_yr, fe_h):
    """The record of a neutron star or black hole of `mass_solar` that formed from a progenitor of `progenitor_solar`: not a star, so
    it has no `params`; its largest mass is the progenitor's, and its radius is its own (its progenitor's giant radius is the history
    rule's, GEN.1012)."""
    return MemberRecord(kind, None, mass_solar, progenitor_solar, radius_au, radius_au, max(progenitor_solar, mass_solar), None,
                        None, age_yr, fe_h, seed)


def of_star(star, seed=0):
    """
    The record of a `Star` that already exists (the primary of a system, which the tree is drawn around): its present values, its
    birth mass as the star holds it, and the largest radius its birth mass and age give. Its `params` are None, as the star is
    built and nothing is built from this record.
    """
    mass = _kg_to_msun(star.mass)
    initial = getattr(star, "initial_mass_sol", None) or mass
    kind = _kind_of(star.yerkes_class)
    radius = _radius_au_of_km(star.radius)
    age_yr = star.age * 1.0e9
    fe_h = star.metallicity_fe_h
    largest = radius if kind in YOUNG_KINDS else _largest_radius_au(initial, age_yr, fe_h, radius)
    return MemberRecord(kind, star.yerkes_class, mass, initial, radius, largest, max(initial, mass), star.temperature,
                        star.luminosity / constants.SOLAR_LUMINOSITY, age_yr, fe_h, seed)
