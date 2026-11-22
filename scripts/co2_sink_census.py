#!/usr/bin/env python3
"""
GEN.1821: what a photosynthesis carbon sink would do to the CO2 of generated worlds (a measurement, never a setting; nothing in `src/` changes).

    python scripts/co2_sink_census.py collect OUT.jsonl SET [--workers 3] [--count N]
    python scripts/co2_sink_census.py report OUT.jsonl

`collect` is `class_census.py collect` with one addition: after `generation.life.apply_life_data` has finished a body, the carbon sink of the body's air
(`CarbonSink.kept_at`) is solved again as if the life it carries buried organic carbon (the report of `reports/co2-carbon-sink-gen1821-2026-11-21.md`):

* the net sink of an oxygenic biosphere is `B = S - W(p)` (`biogenic_oxygen`: the source `S`, the oxygenic weathering `W(p)` at the oxygen of the air; all of `S` while no
  oxygen stands, `R` at the oxygen balance), in moles of carbon a year; a body whose life makes no oxygen has `B = 0` here (classes 2 to 4 have no burial row yet);
* `b` is that over Earth's outgassing per area: `(B / area) / (F_E / area_E)`, with `F_E` = 6.6 Tmol C a year (`waterworld_carbon.SEAFLOOR_SHARE_OF_OUTGASSING`, [Krissansen-Totton and Catling]);
* **option 1** adds `max(b - b_E, 0)` to the excess of the balance (`b_E = R_E / F_E`, Earth's own net sink, 3 / 6.6): life only beyond Earth's;
* **option 2** weights the weathering by `1 - b_E` and adds `b`: the lifeless baseline (an Earth balances with its own biosphere).

The body is then left as it was; the record gets `co2_census`: for the baseline and each option the CO2 pressure (Pa), the surface temperature (K) and the life classes
that still pass their gates, plus `S`, `R`, `p_O2`, `b`.
"""

import argparse
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "src"))

import class_census  # noqa: E402

FLOW_EARTH_MOL_YR = 6.6e12
"""float: Earth's modern outgassing of CO2, mol C a year (the 6.6 Tmol of `waterworld_carbon.SEAFLOOR_SHARE_OF_OUTGASSING`)."""

STATE = {"option": 0, "b": 0.0, "b_earth": 0.0}
"""dict: The option being solved (0 baseline, 1, 2), the body's `b` and Earth's `b_E`, read by the patched balances."""


def extra_term(taken):
    """`(taken, extra)` for the option in force: option 2 weights what the weathering takes by `1 - b_E` and adds `b`; option 1 adds `max(b - b_E, 0)`."""
    if STATE["option"] == 1:
        return taken, max(STATE["b"] - STATE["b_earth"], 0.0)
    if STATE["option"] == 2:
        return taken * (1.0 - STATE["b_earth"]), STATE["b"]
    return taken, 0.0


def install_balances():
    """Replaces the two flux balances of the carbon sink with the same laws plus the extra term (the baseline option 0 changes nothing)."""
    from planetgen.physics import continental_weathering, surface_pressure, waterworld_carbon

    original_land = continental_weathering.pco2_bar
    original_sea = waterworld_carbon.seafloor_pco2_bar

    def land_balance(supply_bar, heat_flux_w_m2, land_fraction, temperature_k, absorbed_flux_w_m2, cycle, seafloor_law=None):
        if STATE["option"] == 0:
            return original_land(supply_bar, heat_flux_w_m2, land_fraction, temperature_k, absorbed_flux_w_m2, cycle, seafloor_law)
        lambda_scale, exponent = cycle
        outgassing = max(heat_flux_w_m2 / surface_pressure.HEAT_FLUX_EARTH_W_M2, 0.0) ** exponent
        if outgassing <= 0.0:
            return original_land(supply_bar, heat_flux_w_m2, land_fraction, temperature_k, absorbed_flux_w_m2, cycle, seafloor_law)
        seafloor_share = waterworld_carbon.SEAFLOOR_SHARE_OF_OUTGASSING if seafloor_law is not None else 0.0
        land = max(land_fraction, 0.0) / surface_pressure.LAND_FRACTION_EARTH

        def excess(p_bar):
            taken = (1.0 - seafloor_share) * land * continental_weathering.rate_over_modern(temperature_k, p_bar, absorbed_flux_w_m2, lambda_scale)
            if seafloor_law is not None:
                taken += seafloor_share * waterworld_carbon.seafloor_rate_over_modern(p_bar, temperature_k, heat_flux_w_m2, seafloor_law)
            taken, extra = extra_term(taken)
            return taken + extra - outgassing

        return surface_pressure.balance_root_bar(supply_bar, excess)

    def sea_balance(supply_bar, heat_flux_w_m2, temperature_of, law, outgassing_exponent):
        if STATE["option"] == 0 or supply_bar <= surface_pressure.MIN_PCO2_BAR:
            return original_sea(supply_bar, heat_flux_w_m2, temperature_of, law, outgassing_exponent)
        outgassing = max(heat_flux_w_m2 / waterworld_carbon.HEAT_FLUX_EARTH_W_M2, 0.0) ** outgassing_exponent
        if outgassing <= 0.0:
            return surface_pressure.MIN_PCO2_BAR

        def excess(p_bar):
            taken = waterworld_carbon.SEAFLOOR_SHARE_OF_OUTGASSING * waterworld_carbon.seafloor_rate_over_modern(p_bar, temperature_of(p_bar), heat_flux_w_m2, law)
            taken, extra = extra_term(taken)
            return taken + extra - outgassing

        return surface_pressure.balance_root_bar(supply_bar, excess)

    continental_weathering.pco2_bar = land_balance
    waterworld_carbon.seafloor_pco2_bar = sea_balance


def biosphere_sink(planet):
    """`(S, R, p_O2, B)` of the oxygenic life of `planet`, mol of carbon a year and Pa; `B = 0` for a body whose life makes no oxygen."""
    from planetgen.physics import atmosphere, biogenic_oxygen, life_air
    living = getattr(getattr(planet, "build", None), "air_living", None) or []
    pairs = life_air._oxygenic(planet, living, forever=False)
    if not pairs or getattr(planet, "body_type", "t") == "g":
        return 0.0, 0.0, 0.0, 0.0
    area = biogenic_oxygen.area_m2(planet)
    source = biogenic_oxygen.source_of(planet, area, [(share, band) for _, share, band in pairs])
    reductant = biogenic_oxygen.reductant_flux(planet, area)
    pascal = atmosphere.gas_partials_pa(planet).get("O2", 0.0)
    weathering = biogenic_oxygen.EARTH_SURPLUS_MOL_YR * (area / biogenic_oxygen.EARTH_AREA_M2) * (pascal / biogenic_oxygen.EARTH_O2_PA) ** biogenic_oxygen.OXIDATIVE_WEATHERING_EXPONENT
    return source, reductant, pascal, max(source - min(weathering, source), 0.0)


def state_of(planet):
    """What the report compares: the CO2 pressure (Pa), the surface temperature (K) and the life classes that pass their gates."""
    from planetgen.physics import atmosphere, life_classes
    scores = life_classes.score_classes(planet)
    return {"co2_pa": atmosphere.gas_partials_pa(planet).get("CO2", 0.0), "t_k": planet.surface_temperature,
            "classes": sorted(class_id for class_id, one in scores.items() if one.score > 0)}


def resolve(planet):
    """Solves the CO2 of `planet` again under each option and leaves the body as it was; returns the record or `None` for a body with no carbon sink."""
    from planetgen.generation.build_record import record_of
    from planetgen.physics import atmosphere, biogenic_oxygen
    sink = getattr(record_of(planet), "co2_sink", None)
    if sink is None or not sink.solved_bar > 0.0:
        return None
    source, reductant, pascal, burial = biosphere_sink(planet)
    area = biogenic_oxygen.area_m2(planet)
    per_area_earth = FLOW_EARTH_MOL_YR / biogenic_oxygen.EARTH_AREA_M2
    STATE["b"] = (burial / area) / per_area_earth
    STATE["b_earth"] = biogenic_oxygen.EARTH_REDUCTANT_MOL_YR / FLOW_EARTH_MOL_YR
    air0 = dict(atmosphere.gas_partials_pa(planet))
    started = (planet.surface_temperature, getattr(planet, "greenhouse_factor", None))
    out = {"S": source, "R": reductant, "p_o2": pascal, "B": burial, "b": STATE["b"], "b_earth": STATE["b_earth"], "base": state_of(planet)}
    from planetgen.physics import biogenic_oxygen as oxygen
    for option in (0, 1, 2):      # option 0 is the re-solve alone with the law as it is: the options are compared with it, since the re-solve itself moves bodies whose stored CO2 was changed after the sink was solved
        STATE["option"] = option
        try:
            changed = oxygen.rebalance_co2(planet)
            out[f"option{option}"] = dict(state_of(planet), changed=bool(changed))
        except Exception as exc:  # noqa: BLE001 -- a body that cannot be solved is a row, not a stop
            out[f"option{option}"] = {"error": f"{type(exc).__name__}: {exc}"[:120]}
        finally:
            STATE["option"] = 0
            planet.surface_temperature, planet.greenhouse_factor = started
            oxygen.store_air(planet, air0)
    return out


def install_hook():
    """After `apply_life_data` of every body, attaches `co2_census` and lets `class_census` record it."""
    from planetgen.generation import life
    original = life.apply_life_data

    def wrapped(planet, keep_stored=False):
        original(planet, keep_stored)
        try:
            planet.co2_census = resolve(planet)
        except Exception as exc:  # noqa: BLE001
            planet.co2_census = {"error": f"{type(exc).__name__}: {exc}"[:120]}

    life.apply_life_data = wrapped
    class_census.BODY_ATTRS.append("co2_census")


def collect(args):
    install_balances()
    install_hook()
    class_census.CLASS_M = None
    class_census.draw_sample(args.out, args.set, args.workers, args.count)


def report(args):
    import collections
    import gzip
    import json
    import statistics
    rows = []
    opener = gzip.open if args.path.endswith(".gz") else open
    with opener(args.path, "rt") as handle:
        for line in handle:
            try:
                record = json.loads(line)
            except ValueError:
                continue
            for body in record.get("bodies", []):
                census = body.get("co2_census")
                if isinstance(census, dict) and "base" in census:
                    rows.append((body, census))
    print(f"{len(rows)} bodies with a carbon sink\n")
    groups = {"no life": lambda b, c: b.get("life_class") is None,
              "life, no oxygen standing": lambda b, c: b.get("life_class") is not None and c["p_o2"] < 100.0,
              "life, oxygen risen": lambda b, c: b.get("life_class") is not None and c["p_o2"] >= 100.0}
    for option in (0, 1, 2):
        print(f"## {'Re-solve alone, compared with the stored air' if option == 0 else 'Option ' + str(option) + ', compared with the re-solve alone'}\n\n| group | bodies | CO2 x2 or more (up) | CO2 half or less (down) | median dT K | p90 |dT| K | cross 273 K | cross 373 K | CO2 under 0.5 Pa | class set changed |\n|---|---|---|---|---|---|---|---|---|---|---|")
        for name, test in groups.items():
            members = [(b, c) for b, c in rows if test(b, c) and isinstance(c.get(f"option{option}"), dict) and "co2_pa" in c[f"option{option}"]]
            if not members:
                print(f"| {name} | 0 | | | | | | | | | |")
                continue
            up = down = freeze = boil = floor = changed = 0
            deltas = []
            for body, census in members:
                before, after = (census["base"] if option == 0 else census.get("option0", census["base"])), census[f"option{option}"]
                if before["co2_pa"] > 0:
                    ratio = after["co2_pa"] / before["co2_pa"]
                    up += ratio >= 2.0
                    down += ratio <= 0.5
                delta = after["t_k"] - before["t_k"]
                deltas.append(delta)
                freeze += (before["t_k"] >= 273.15) != (after["t_k"] >= 273.15)
                boil += (before["t_k"] >= 373.15) != (after["t_k"] >= 373.15)
                floor += after["co2_pa"] < 0.5 <= before["co2_pa"]
                changed += before["classes"] != after["classes"]
            deltas.sort()
            p90 = sorted(abs(d) for d in deltas)[int(0.9 * (len(deltas) - 1))]
            print(f"| {name} | {len(members)} | {up} | {down} | {statistics.median(deltas):.1f} | {p90:.1f} | {freeze} | {boil} | {floor} | {changed} |")
        print()


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    commands = parser.add_subparsers(dest="command", required=True)
    drawer = commands.add_parser("collect")
    drawer.add_argument("out")
    drawer.add_argument("set", choices=class_census.SETS)
    drawer.add_argument("--workers", type=int, default=3)
    drawer.add_argument("--count", type=int)
    reporter = commands.add_parser("report")
    reporter.add_argument("path")
    args = parser.parse_args()
    if args.command == "collect":
        collect(args)
    else:
        report(args)


if __name__ == "__main__":
    main()
