# tests/test_star_composition_sweep.py

"""
TEST.211: a star composition sweep. Stars of every kind (protostar to black hole, brown dwarfs and remnants included) are made the way
the generator makes them, and one set of rules is run over their compositions (`planetgen.physics.stellar_composition`):

- the fractions add up to one and stay in range (design note part A, invariant 1, and the helium line `Y = Y_P + g Z`, invariant 3);
- they follow [Fe/H] and age as the relations say (metals rise with iron, the high-alpha sequence waits for the age gate and
  takes the metal-poor old stars);
- a red giant has its first dredge-up and no other kind changes its surface;
- remnants have their layers (a white dwarf its core, atmosphere, envelope and radius; a neutron star its crust; a black hole its horizon),
  a Wolf-Rayet star is hydrogen-free, a brown dwarf reports its burning limits;
- no star lacks a composition;
- the result is the same in another order, again, and at four workers as at one.

The stars come from two sources: seeded sectors filled on the census harness (`tests.audit.census`, the rates of the rare kinds
boosted as the representation audit boosts them), which are the generated galaxy as the generator fills it, and a ladder of masses
and ages through the generator's own stellar evolution (`evolve_star`, `wolf_rayet.state`) that makes sure every kind is present.
The rules are plain functions of a star and its composition, so `test_each_rule_fails_when_it_is_broken` breaks them one at a time
on copies of the compositions and checks that the rule named for it, and no other, reports it. The sweep reads the downloaded
data (the isochrones, the white-dwarf sequences, the neutron-star crust), like the other tests of the generator's own stars.
"""

import concurrent.futures
import dataclasses
import math
import multiprocessing
import random

import pytest

from planetgen import tuning
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
from planetgen.physics import brown_dwarf, constants, wolf_rayet
from planetgen.physics import stellar_abundances as sa
from planetgen.physics import stellar_composition as sc
from planetgen.physics import post_main_sequence
from planetgen.physics.stellar_evolution import evolve_star, main_sequence_lifetime_gy
from planetgen.util import draw
from tests.audit import census as audit
from tests.audit import representation as rep

SECTOR_REGIONS = ("thin-disk-arm", "bulge", "halo", "globular-cluster", "satellite")
"""tuple: The audit regions whose seeded sectors the sweep fills."""
SYSTEMS_PER_SECTOR = 32
LADDER_MASSES = (0.04, 0.06, 0.09, 0.3, 0.8, 1.0, 1.5, 2.0, 2.5, 5.0, 9.0, 15.0, 25.0, 40.0, 80.0)
LADDER_AGES = (0.0002, 0.005, 0.05, 0.5, 2.0, 5.0, 9.0, 11.0, 13.0)
LADDER_LIFETIME_SHARES = (0.3, 0.95, 1.0, 1.03, 1.06, 1.1, 1.2, 1.4, 2.0)
AGB_MASSES = (1.0, 1.5, 2.0, 2.5)
AGB_FE_H = (-1.0, -0.2, 0.2)
SOLAR_DISTANCE_LY = 26000.0
OLD = tuning.STELLAR_THICK_START_GY
"""float: The age gate of the high-alpha sequence, Gy."""


@dataclasses.dataclass(frozen=True)
class Star:
    """What a star's composition is derived from, and where it came from."""

    source: str
    kind: str
    fe_h: float | None
    age_gy: float
    mass_sol: float
    key: str
    subtype: str | None = None
    teff_k: float | None = None
    spin: float = 0.0


def makeup_of(star):
    """The composition of `star`."""
    return sc.star_makeup(star.kind, star.fe_h, star.age_gy, star.mass_sol, star.key, wolf_rayet_subtype=star.subtype,
                          spin=star.spin, teff_k=star.teff_k)


# --- The stars -------------------------------------------------------------------------------------------------------------------

def _from_sector_star(star, index, region):
    """A star of a filled sector as the generator stored it."""
    if isinstance(star, NeutronStar):
        return Star("sector", sc.NEUTRON_STAR, None, star.age, getattr(star, "mass_solar", 1.4), f"{region}:{index}")
    if isinstance(star, BlackHole):
        return Star("sector", sc.BLACK_HOLE, None, star.age, getattr(star, "mass_solar", 10.0), f"{region}:{index}",
                    spin=min(max(getattr(star, "spin", 0.0) or 0.0, 0.0), 0.99))
    mass_sol = star.mass / constants.SOLAR_MASS_TO_KG
    kind = sc.kind_of(star.yerkes_class, mass_sol, None, sa.PROTO_Z * 10.0 ** star.metallicity_fe_h)
    return Star("sector", kind, star.metallicity_fe_h, star.age, mass_sol,
                sc.abundance_key(mass=star.mass, distance_ly=star.galactic_center_dist_ly),
                teff_k=star.temperature)


def sector_stars(patch):
    """Every star of a few seeded sectors in each region class, filled with the rare kinds' rates boosted."""
    rep.boosted_rates(patch)
    found = []
    for name in SECTOR_REGIONS:
        case = dataclasses.replace(rep.region_cases((name,), 1)[0], num_systems=SYSTEMS_PER_SECTOR)
        sector = audit.fill_region_sector(case, audit.seed_for(name, 0))
        for entry in sector.entries:
            for number, star in enumerate(entry.star_system.stars):
                found.append(_from_sector_star(star, f"{len(found)}:{number}", name))
    return found


def middle_of_the_asymptotic_giant_branch_gy(mass, fe_h):
    """The age half-way through the asymptotic-giant phase of a star of `mass` and `fe_h` (a window a few hundredths of a Gy wide, after
    the horizontal branch), from the same times `evolve_star` uses."""
    window = post_main_sequence.timeline(mass, fe_h).window(post_main_sequence.ASYMPTOTIC_GIANT)
    return 0.5 * (window.start_gy + window.end_gy)


def ladder_stars():
    """Stars of masses from 0.04 to 80 solar masses at ages from before the main sequence to long after it, through the generator's
    stellar evolution; the Wolf-Rayet stars from its own Wolf-Rayet state."""
    stream = draw.Stream("test211-ladder")
    found = []
    for mass in LADDER_MASSES:
        lifetime = main_sequence_lifetime_gy(mass)
        for age in LADDER_AGES + tuple(lifetime * share for share in LADDER_LIFETIME_SHARES):
            fe_h = -2.0 + 2.4 * stream.random()
            key = sc.abundance_key(mass=mass, distance_ly=SOLAR_DISTANCE_LY)
            if mass < sa.hydrogen_burning_limit_for_fe_h(fe_h):
                # GEN.936: under the burning limit of its place an object has no main sequence; it is a brown dwarf on its cooling track
                found.append(Star("ladder", sc.BROWN_DWARF, fe_h, age, mass, key, teff_k=brown_dwarf.state(mass, age)["temperature_k"]))
                continue
            state = evolve_star(mass, age, rng=stream, fe_h=fe_h)
            if state is None:
                if mass < 25.0:
                    found.append(Star("ladder", sc.NEUTRON_STAR, None, age, 1.4, key))
                else:
                    found.append(Star("ladder", sc.BLACK_HOLE, None, age, mass, key, spin=stream.random() * 0.99))
                continue
            kind = sc.kind_of(state["yerkes_class"], state["mass_sol"], state.get("phase"))
            found.append(Star("ladder", kind, fe_h, age, state["mass_sol"], key, subtype=state.get("wr_subtype"),
                              teff_k=state["temperature_k"]))
    for mass in AGB_MASSES:
        for fe_h in AGB_FE_H:
            age = middle_of_the_asymptotic_giant_branch_gy(mass, fe_h)
            state = evolve_star(mass, age, rng=stream, fe_h=fe_h)
            key = sc.abundance_key(mass=mass, distance_ly=SOLAR_DISTANCE_LY)
            found.append(Star("ladder", sc.kind_of(state["yerkes_class"], state["mass_sol"], state.get("phase")), fe_h, age,
                              state["mass_sol"], key, teff_k=state["temperature_k"]))
    for index in range(60):
        age = 0.003 + 0.001 * stream.random()
        state = wolf_rayet.state(40.0 + 40.0 * stream.random(), 15.0 + 20.0 * stream.random(), age, 0.003, 0.004, stream)
        found.append(Star("ladder", sc.WOLF_RAYET, -2.0 + 2.4 * stream.random(), age, state["mass_sol"],
                          sc.abundance_key(mass=state["mass_sol"], distance_ly=SOLAR_DISTANCE_LY + index),
                          subtype=state["wr_subtype"], teff_k=state["temperature_k"]))
    return found


@pytest.fixture(scope="module")
def galaxy():
    """`(stars, compositions)`: every star of the sweep and its composition."""
    patch = pytest.MonkeyPatch()
    try:
        stars = sector_stars(patch) + ladder_stars()
    finally:
        patch.undo()
    return stars, [makeup_of(star) for star in stars]


# --- The rules -------------------------------------------------------------------------------------------------------------------
# A star rule takes a star and its composition and returns a sentence for what is wrong, or None. Names are the rules' own.

def _compositions(makeup):
    return [composition for composition in (makeup.birth, makeup.surface) if composition is not None]


def rule_fractions_add_up(star, makeup):
    for composition in _compositions(makeup):
        if abs(math.fsum(composition.elements.values()) - 1.0) > 1e-9 or abs(composition.x + composition.y + composition.z - 1.0) > 1e-9:
            return "the element fractions or X + Y + Z do not add up to one"
    fractions = makeup.detail.get("surface_fractions")
    if fractions is not None and abs(math.fsum(fractions.values()) - 1.0) > 1e-9:
        return "the Wolf-Rayet surface fractions do not add up to one"
    layers = makeup.detail.get("layers")
    if layers is not None and abs(math.fsum(layers.values()) - 1.0) > 1e-3:
        return "the white dwarf's layers do not add up to one"
    return None


def rule_fractions_in_range(star, makeup):
    for composition in _compositions(makeup):
        if min(composition.elements.values()) < 0.0 or max(composition.elements.values()) > 1.0:
            return "an element fraction is outside 0 to 1"
        if not (0.5 < composition.x < 0.8 and sa.PRIMORDIAL_Y - 1e-12 <= composition.y < 0.5 and 0.0 <= composition.z < 0.1):
            return "X, Y or Z is outside the range a star of this [Fe/H] can have"
    fractions = makeup.detail.get("surface_fractions")
    if fractions is not None and (min(fractions.values()) < 0.0 or max(fractions.values()) > 1.0):
        return "a Wolf-Rayet surface fraction is outside 0 to 1"
    return None


def rule_helium_follows_the_metals(star, makeup):
    """Part A invariant 3: the birth helium is the primordial helium plus `g` times the metals."""
    if makeup.birth is not None and abs(makeup.birth.y - sa.PRIMORDIAL_Y - sa.HELIUM_PER_METAL * makeup.birth.z) > 1e-9:
        return "the birth helium is not Y_P + g Z"
    return None


def rule_the_branch_waits_for_the_age_gate(star, makeup):
    if makeup.kind in sc.WITH_BIRTH_COMPOSITION:
        if makeup.branch not in sa.BRANCHES:
            return "no abundance sequence"
        if star.age_gy < OLD and makeup.branch != sa.BRANCH_THIN:
            return "a star younger than the age gate is on the high-alpha sequence"
    elif makeup.branch is not None:
        return "a remnant with no metallicity has an abundance sequence"
    return None


def rule_only_a_red_giant_has_a_dredge_up(star, makeup):
    if makeup.birth is None:
        return None
    if makeup.kind == sc.RED_GIANT:
        birth, surface = makeup.birth, makeup.surface
        if not (surface.elements["N"] > birth.elements["N"] and surface.elements["C"] < birth.elements["C"] and surface.y > birth.y):
            return "a red giant's surface has no first dredge-up (more nitrogen, less carbon, more helium)"
    elif makeup.surface != makeup.birth:
        return f"a {makeup.kind} star's surface differs from its birth composition"
    return None


def rule_every_star_has_a_composition(star, makeup):
    if makeup.kind != star.kind:
        return "the composition is of another kind"
    if makeup.kind in sc.WITH_BIRTH_COMPOSITION and (makeup.birth is None or makeup.surface is None or makeup.birth.fe_h != star.fe_h):
        return "a star with a metallicity has no birth composition at it"
    if makeup.kind == sc.NEUTRON_STAR and not makeup.detail.get("outer_crust"):
        return "a neutron star has no outer crust"
    if makeup.kind == sc.BLACK_HOLE and makeup.detail.get("hole") is None:
        return "a black hole has no horizon quantities"
    if makeup.kind == sc.WOLF_RAYET and not makeup.detail.get("surface_fractions"):
        return "a Wolf-Rayet star has no surface fractions"
    if makeup.kind == sc.BROWN_DWARF and "burns_hydrogen" not in makeup.detail:
        return "a brown dwarf reports no burning limits"
    return None


def rule_remnants_have_their_layers(star, makeup):
    detail = makeup.detail
    if makeup.kind == sc.WHITE_DWARF:
        if detail.get("core") not in ("helium", "carbon-oxygen", "oxygen-neon") or detail.get("atmosphere") not in ("DA", "non-DA"):
            return "a white dwarf has no core or atmosphere"
        if detail.get("envelope") != ("thick" if detail["atmosphere"] == "DA" else "thin"):
            return "a white dwarf's envelope does not follow its atmosphere"
        if not detail.get("radius_solar") or not 0.003 < detail["radius_solar"] < 0.05:
            return "a white dwarf has no radius, or one no white dwarf has"
        if detail["core"] == "carbon-oxygen" and detail.get("layers") is None:
            return "a carbon-oxygen white dwarf has no layers"
        if detail["polluted"] and star.teff_k >= constants.COOL_WHITE_DWARF_BELOW_K:
            return "a hot white dwarf is polluted"
    elif makeup.kind == sc.NEUTRON_STAR:
        if not 2.0e11 < detail["neutron_drip_g_cm3"] < 8.0e11:
            return "a neutron star's crust does not end at neutron drip"
        if not all(isinstance(nucleus, str) and "-" in nucleus for nucleus in detail["outer_crust"]):
            return "a neutron star's crust is not a list of nuclei"
    elif makeup.kind == sc.BLACK_HOLE:
        hole = detail["hole"]
        expected_km = 2.0 * constants.G * star.mass_sol * constants.SOLAR_MASS_TO_KG / constants.SPEED_OF_LIGHT_M_S ** 2 / 1e3
        if abs(hole.schwarzschild_radius_km - expected_km) > 1e-6 * expected_km:
            return "a black hole's Schwarzschild radius is not 2GM/c^2"
    elif makeup.kind == sc.WOLF_RAYET:
        if detail["surface_fractions"].get("H", 0.0) > 0.0 and sc.wolf_rayet_row(star.subtype, star.teff_k) != "WN-late":
            return "a Wolf-Rayet star that is not a late WN star has hydrogen"
    elif makeup.kind == sc.BROWN_DWARF:
        if detail["burns_hydrogen"] != (star.mass_sol > sa.hydrogen_burning_limit_sol(makeup.birth.z)):
            return "a brown dwarf's hydrogen burning does not follow its mass"
        if detail["burned_deuterium"] != (detail["mass_jupiter"] > detail["deuterium_limit_mj"]):
            return "a brown dwarf's deuterium burning does not follow its mass"
    return None


STAR_RULES = {
    "fractions add up": rule_fractions_add_up,
    "fractions in range": rule_fractions_in_range,
    "helium follows the metals": rule_helium_follows_the_metals,
    "the branch waits for the age gate": rule_the_branch_waits_for_the_age_gate,
    "only a red giant has a dredge-up": rule_only_a_red_giant_has_a_dredge_up,
    "every star has a composition": rule_every_star_has_a_composition,
    "remnants have their layers": rule_remnants_have_their_layers,
}
"""dict: The rules of one star, by name."""


def rank_correlation(xs, ys):
    """Spearman's rank correlation of two lists (no ties handled: they are floats)."""
    def ranks(values):
        order = sorted(range(len(values)), key=values.__getitem__)
        out = [0.0] * len(values)
        for rank, index in enumerate(order):
            out[index] = float(rank)
        return out
    rx, ry = ranks(xs), ranks(ys)
    mean = (len(xs) - 1) / 2.0
    covariance = sum((a - mean) * (b - mean) for a, b in zip(rx, ry))
    return covariance / math.sqrt(sum((a - mean) ** 2 for a in rx) * sum((b - mean) ** 2 for b in ry))


def rule_the_metals_follow_the_iron(stars, makeups):
    """Over the stars with a metallicity, the metal fraction rises with [Fe/H]."""
    pairs = [(star.fe_h, makeup.birth.z) for star, makeup in zip(stars, makeups) if makeup.birth is not None]
    if len(pairs) < 50:
        return "too few stars with a metallicity to see the trend"
    if rank_correlation([p[0] for p in pairs], [p[1] for p in pairs]) < 0.9:
        return "the metal fraction does not follow [Fe/H]"
    return None


def rule_the_old_and_poor_are_high_alpha(stars, makeups):
    """Past the age gate the high-alpha sequence takes the metal-poor stars (probability 0.97 at -0.8) and leaves the metal-rich ones (under
    0.01 above 0.2)."""
    poor = [m.branch == sa.BRANCH_HIGH_ALPHA for s, m in zip(stars, makeups) if m.birth is not None and s.age_gy >= OLD and s.fe_h <= -0.8]
    rich = [m.branch == sa.BRANCH_HIGH_ALPHA for s, m in zip(stars, makeups) if m.birth is not None and s.age_gy >= OLD and s.fe_h >= 0.2]
    if len(poor) < 10 or len(rich) < 3:
        return "too few old stars to see the high-alpha sequence"
    if sum(poor) / len(poor) < 0.8:
        return "the old metal-poor stars are not mostly high-alpha"
    if sum(rich) / len(rich) > 0.2:
        return "the old metal-rich stars are high-alpha"
    return None


def rule_every_kind_is_present(stars, makeups):
    missing = sorted(set(sc.KINDS) - {star.kind for star in stars})
    return f"the sweep never made {missing}" if missing else None


def rule_the_same_in_another_order_and_again(stars, makeups):
    order = list(range(len(stars)))
    random.Random(211).shuffle(order)
    again = {index: makeup_of(stars[index]) for index in order}
    if any(again[index] != makeups[index] for index in order):
        return "a composition changed with the order the stars were asked for in"
    return None


def _makeup_in_worker(star):
    return makeup_of(star)


def rule_the_same_at_four_workers(stars, makeups):
    # Forked, so the workers keep the cache folder of the data the parent read.
    with concurrent.futures.ProcessPoolExecutor(max_workers=4, mp_context=multiprocessing.get_context("fork")) as pool:
        parallel = list(pool.map(_makeup_in_worker, stars, chunksize=16))
    return None if parallel == makeups else "a composition differs between one worker and four"


SET_RULES = {
    "the metals follow the iron": rule_the_metals_follow_the_iron,
    "the old and poor are high-alpha": rule_the_old_and_poor_are_high_alpha,
    "every kind is present": rule_every_kind_is_present,
    "the same in another order and again": rule_the_same_in_another_order_and_again,
    "the same at four workers": rule_the_same_at_four_workers,
}
"""dict: The rules over the whole sweep, by name."""


def problems(stars, makeups, set_rules=SET_RULES):
    """`{rule name: [what is wrong]}` over the sweep: every star rule on every star, then the rules over all of it."""
    found = {}
    for star, makeup in zip(stars, makeups):
        for name, rule in STAR_RULES.items():
            message = rule(star, makeup)
            if message is not None:
                found.setdefault(name, []).append(f"{star.source} {star.kind} [Fe/H] {star.fe_h} age {star.age_gy:.4g}: {message}")
    for name, rule in set_rules.items():
        message = rule(stars, makeups)
        if message is not None:
            found.setdefault(name, []).append(message)
    return found


# --- The sweep -------------------------------------------------------------------------------------------------------------------

def test_the_sweep_has_every_kind_from_a_filled_galaxy_and_from_the_ladder(galaxy):
    stars, _makeups = galaxy
    assert {star.kind for star in stars} == set(sc.KINDS)
    sector = [star for star in stars if star.source == "sector"]
    # Giants are well under 1 percent of the population model's stars (GEN.927, GEN.935), so 78 stars rarely hold one; the ladder gives every kind.
    # The filled sectors must still hold more than one kind (main-sequence stars and remnants).
    assert len(sector) >= 20 and len({star.kind for star in sector}) >= 2, "the filled sectors hold a spread of stars"
    assert all(star.source in ("sector", "ladder") for star in stars)


def test_no_rule_finds_anything_wrong_with_the_generators_stars(galaxy):
    stars, makeups = galaxy
    found = problems(stars, makeups)
    assert found == {}, "\n".join(f"{name}: {len(messages)}, e.g. {messages[0]}" for name, messages in found.items())


def test_every_star_of_a_filled_sector_has_a_composition(galaxy):
    stars, makeups = galaxy
    for star, makeup in zip(stars, makeups):
        if star.source == "sector":
            assert makeup.kind == star.kind and (makeup.birth is not None or makeup.kind in (sc.NEUTRON_STAR, sc.BLACK_HOLE)), star


def test_the_metals_rise_with_the_iron_star_by_star_on_the_same_draws(galaxy):
    """The relation, not just the trend: the same star at a higher [Fe/H] with the same key and the same sequence has more metal."""
    stars, _makeups = galaxy
    checked = 0
    for star in stars:
        if star.kind not in sc.WITH_BIRTH_COMPOSITION or star.source != "ladder":
            continue
        low = makeup_of(star)
        high = makeup_of(dataclasses.replace(star, fe_h=star.fe_h + 0.3))
        if low.branch == high.branch:
            assert high.birth.z > low.birth.z, star
            checked += 1
    assert checked > 100


# --- The rules fail when they are broken -----------------------------------------------------------------------------------------

def _tamper(makeups, picks, change):
    """A copy of `makeups` with `change(makeup)` applied to the ones `picks(index, makeup)` selects."""
    return [change(makeup) if picks(index, makeup) else makeup for index, makeup in enumerate(makeups)]


def _with_elements(composition, **changes):
    return dataclasses.replace(composition, elements={**composition.elements, **changes})


def _first(makeups, kind):
    return next(index for index, makeup in enumerate(makeups) if makeup.kind == kind)


def _break_the_sums(makeup):
    birth = _with_elements(makeup.birth, H=makeup.birth.elements["H"] + 0.01)
    return dataclasses.replace(makeup, birth=birth, surface=birth)


def _break_the_range(makeup):
    birth = _with_elements(makeup.birth, H=makeup.birth.elements["H"] + 0.001, Fe=-0.001)
    return dataclasses.replace(makeup, birth=birth, surface=birth)


def _break_the_helium(makeup):
    birth = dataclasses.replace(makeup.birth, x=makeup.birth.x - 0.01, y=makeup.birth.y + 0.01)
    return dataclasses.replace(makeup, birth=birth, surface=birth)


def _break_the_gate(makeup):
    return dataclasses.replace(makeup, branch=sa.BRANCH_HIGH_ALPHA)


def _break_the_dredge_up(makeup):
    return dataclasses.replace(makeup, surface=makeup.birth)


def _add_a_dredge_up(makeup):
    return dataclasses.replace(makeup, surface=_with_elements(makeup.birth, N=makeup.birth.elements["N"] * 2.0))


def _lose_the_composition(makeup):
    return dataclasses.replace(makeup, birth=None)


def _lose_the_crust(makeup):
    return dataclasses.replace(makeup, detail={**makeup.detail, "outer_crust": ()})


def _lose_the_atmosphere(makeup):
    return dataclasses.replace(makeup, detail={**makeup.detail, "atmosphere": None})


def _warp_the_horizon(makeup):
    hole = dataclasses.replace(makeup.detail["hole"], schwarzschild_radius_km=makeup.detail["hole"].schwarzschild_radius_km * 1.01)
    return dataclasses.replace(makeup, detail={**makeup.detail, "hole": hole})


def _give_it_hydrogen(makeup):
    fractions = {**makeup.detail["surface_fractions"], "H": 0.1}
    return dataclasses.replace(makeup, detail={**makeup.detail, "surface_fractions": fractions})


def _flip_the_burning(makeup):
    return dataclasses.replace(makeup, detail={**makeup.detail, "burns_hydrogen": not makeup.detail["burns_hydrogen"]})


@pytest.mark.parametrize("name, kind, change", [
    ("fractions add up", sc.MAIN_SEQUENCE, _break_the_sums),
    ("fractions in range", sc.MAIN_SEQUENCE, _break_the_range),
    ("helium follows the metals", sc.MAIN_SEQUENCE, _break_the_helium),
    ("the branch waits for the age gate", sc.PRE_MAIN_SEQUENCE, _break_the_gate),
    ("only a red giant has a dredge-up", sc.RED_GIANT, _break_the_dredge_up),
    ("only a red giant has a dredge-up", sc.SUBGIANT, _add_a_dredge_up),
    ("every star has a composition", sc.MAIN_SEQUENCE, _lose_the_composition),
    ("every star has a composition", sc.NEUTRON_STAR, _lose_the_crust),
    ("remnants have their layers", sc.WHITE_DWARF, _lose_the_atmosphere),
    ("remnants have their layers", sc.BLACK_HOLE, _warp_the_horizon),
    ("remnants have their layers", sc.WOLF_RAYET, _give_it_hydrogen),
    ("remnants have their layers", sc.BROWN_DWARF, _flip_the_burning),
])
def test_each_star_rule_fails_when_it_is_broken_and_only_that_rule(galaxy, name, kind, change):
    stars, makeups = galaxy
    index = next(index for index, (star, makeup) in enumerate(zip(stars, makeups))
                 if makeup.kind == kind and (name != "the branch waits for the age gate" or star.age_gy < OLD)
                 and (kind != sc.WOLF_RAYET or sc.wolf_rayet_row(star.subtype, star.teff_k) != "WN-late"))
    broken = list(makeups)
    broken[index] = change(makeups[index])
    found = problems(stars, broken, set_rules={})
    assert name in found, f"{name} did not see the break"
    assert all(message.startswith(f"{stars[index].source} {kind}") for message in found[name])
    assert set(found) <= {name, "fractions add up", "fractions in range", "helium follows the metals"}, found.keys()
    assert problems(stars, makeups, set_rules={}) == {}


def test_the_set_rules_fail_when_they_are_broken(galaxy):
    stars, makeups = galaxy
    flat = [dataclasses.replace(m, birth=dataclasses.replace(m.birth, z=0.0)) if m.birth is not None else m for m in makeups]
    assert "the metals follow the iron" in problems(stars, flat)
    always_thin = [dataclasses.replace(m, branch=sa.BRANCH_THIN) if m.birth is not None else m for m in makeups]
    assert "the old and poor are high-alpha" in problems(stars, always_thin)
    assert "every kind is present" in problems([s for s in stars if s.kind != sc.WOLF_RAYET], [m for m in makeups if m.kind != sc.WOLF_RAYET])
    swapped = list(makeups)
    swapped[0], swapped[1] = swapped[1], swapped[0]
    found = problems(stars, swapped)
    assert "the same in another order and again" in found and "the same at four workers" in found


# --- The distribution ------------------------------------------------------------------------------------------------------------

def test_the_solar_neighbourhood_mix_of_stars_has_the_ratios_the_survey_measured():
    """A-V2: the stars of a solar-neighbourhood-like draw ([Fe/H] normal about -0.10 with 0.22 dex, ages to 12 Gy), taken through the
    generator's own derivation, against the medians of Brewer and Fischer 2016 (Fe/Si 1.81 within 0.2, Mg/Si 0.93 within 0.1, C/O 0.385
    within 0.05) and inside the Thiabaud 2015 box (Fe/Si 1.34 to 2.81, Mg/Si 0.66 to 1.85, C/O 0.12 to 0.75) for 80 percent or more."""
    generator = random.Random(2026)
    ratios = []
    for index in range(20000):
        fe_h = min(max(generator.gauss(-0.10, 0.22), -3.5), 0.6)
        makeup = sc.star_makeup(sc.SUBGIANT, fe_h, 12.0 * generator.random(), 1.0, f"test211-dist:{index}")
        ratios.append((sa.mass_ratio(makeup.birth, "Fe", "Si"), sa.mass_ratio(makeup.birth, "Mg", "Si"), sa.mass_ratio(makeup.birth, "C", "O")))
    median = lambda values: sorted(values)[len(values) // 2]  # noqa: E731
    assert median([r[0] for r in ratios]) == pytest.approx(1.81, abs=0.2)
    assert median([r[1] for r in ratios]) == pytest.approx(0.93, abs=0.1)
    assert median([r[2] for r in ratios]) == pytest.approx(0.385, abs=0.05)
    inside = sum(1.34 <= fe <= 2.81 and 0.66 <= mg <= 1.85 and 0.12 <= c <= 0.75 for fe, mg, c in ratios) / len(ratios)
    assert inside >= 0.8
