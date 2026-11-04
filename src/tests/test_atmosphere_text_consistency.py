"""
GEN.666: the atmosphere phrase and the component table of a page name the same gases.

A class B planet's page said "a mix of helium, sodium, and oxygen" above a table of "Oxygen 100%": the helium and sodium of the
class's balance were named in the phrase and held in the surface pressure, but only the ten stored gases were in the table.
Since GEN.424 and GEN.411 every gas of the air is stored (there is no class balance), and the table holds them all, as shares of the
whole surface pressure.

GEN.415: the text is written from the stored mixture alone, so every gas it names is one the body holds above the stated share, under
the gas table's name; the last test builds a planet of every class and zone.
"""

import re

import pytest

from planetgen import tuning
from planetgen.describe import build as describe_build
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import atmosphere as atmosphere_physics
from planetgen.physics import gases

from tests.fuzz_support import deterministic_entropy
from tests.test_planets import ZONE_SEED_IDS, ZONE_SEEDS, host_star, make_planet  # noqa: F401 -- the shared G2V host and the planet builder


def _bodies(seeds):
    for seed in seeds:
        with deterministic_entropy(9000 + seed):
            system = StarSystem(system_config=SystemConfig())
        for planet in system.planets + system.secondary_planets:
            if planet.body_type == "a":
                continue
            yield planet
            yield from getattr(planet, "moons", ())


def _table(body):
    block = describe_build._atmosphere_composition(body)
    return {part["name"]: part["share"] for part in block["parts"]}


def _phrase_gases(body):
    names = [atmosphere_physics.species_name(gas) for gas in atmosphere_physics.SPECIES] + ["helium", "sodium", "ammonia"]
    names += [atmosphere_physics.gas_name(formula) for formula in ("HCl", "HF", "N2O", "PH3", "NO2", "CH3Cl")]       # the trace gases of GEN.524
    text, found = body.atmosphere, set()
    for name in sorted(names, key=len, reverse=True):  # "hydrogen sulfide" before "hydrogen"
        if name in text:
            found.add(name)
            text = text.replace(name, "")
    return found


def test_the_table_has_every_gas_the_phrase_names_and_its_shares_make_one():
    seen = 0
    for body in _bodies(range(80)):
        if body.body_type == "g" or body.atmosphere in (None, "None"):
            continue
        table = _table(body)
        assert abs(sum(table.values()) - 1.0) < 1e-9, (body.name, table)
        phrase = _phrase_gases(body)
        assert phrase <= {name.lower() for name in table}, (body.name, body.atmosphere, table)
        seen += 1
    assert seen >= 15       # about a quarter of the bodies of 80 systems hold air since GEN.411 solves it from the body (it was 21 before)


@pytest.mark.parametrize("zone,seed", ZONE_SEEDS, ids=ZONE_SEED_IDS)
def test_every_planet_says_only_the_gases_it_holds_above_the_stated_share(host_star, zone, seed):
    """GEN.415: for a planet of every zone the text lists the stored gases from `MAJOR_FRACTION` (the four largest) and the
    traces from `TRACE_FRACTION` (the two largest), names each as the gas table does, and no other; the composition table of the
    page holds the same gases and adds up to the surface pressure."""
    for extra in range(3):
        planet = make_planet(host_star, zone, seed * 3 + extra)
        cls = planet.planet_class
        text = planet.atmosphere
        assert text != atmosphere_physics.PENDING
        if planet.atmospheric_pressure == 0.0:
            assert text == "None"
            continue
        air = atmosphere_physics.gas_partials_pa(planet)
        total = sum(air.values())
        assert total == pytest.approx(planet.atmospheric_pressure, rel=1e-9)
        shares = {gases.gas(formula).name: value / total for formula, value in air.items()}
        named = _text_names(text)
        assert named and set(named) <= set(shares), (cls, zone, text, shares)
        majors = sorted((name for name, share in shares.items() if share >= atmosphere_physics.MAJOR_FRACTION), key=lambda n: -shares[n])[:4]
        traces = sorted((name for name, share in shares.items()
                         if atmosphere_physics.TRACE_FRACTION <= share < atmosphere_physics.MAJOR_FRACTION), key=lambda n: -shares[n])[:2]
        assert named == majors + traces, (cls, zone, text, shares)
        table = _table(planet)
        assert {name.lower() for name in table} == set(shares), (cls, zone, table, shares)


def _text_names(text):
    """The gas names an atmosphere text lists, in order: the mix and then the traces, without the adjectives."""
    mix, _, traces = re.sub(r"^a .*?mix of (mostly )?", "", text).partition(", with traces of ")
    return [name for part in (mix, traces) for name in re.split(r", and |, | and ", part) if name]
