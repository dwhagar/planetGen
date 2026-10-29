#!/usr/bin/env python3
"""
Builds the readable table of the life class rules (GEN.886) from the rule file the classifier reads.

    python scripts/build_life_class_rules.py            # write docs/design/habitability/life-class-rules.md
    python scripts/build_life_class_rules.py --check    # build in memory and fail if the committed file differs

The rules are `src/planetgen/physics/data/life_classes.json`; this page is copied from it and never typed by hand, so it cannot drift.
`src/tests/test_life_class_rules_doc.py` fails when it does.
"""

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RULES = os.path.join(HERE, "..", "src", "planetgen", "physics", "data", "life_classes.json")
OUT = os.path.join(HERE, "..", "docs", "design", "habitability", "life-class-rules.md")

GATE_WORDS = (("min", ">="), ("max", "<="), ("gt", ">"), ("lt", "<"), ("equals", "="), ("not_equals", "!="), ("in", "in"), ("not_in", "not in"))


def _cell(text):
    return str(text).replace("|", "\\|").replace("\n", " ")


def gate_text(gate):
    for key, word in GATE_WORDS:
        if key in gate:
            value = gate[key]
            if isinstance(value, list):
                value = ", ".join(str(item) for item in value)
            return f"{gate['quantity']} {word} {value}"
    return gate["quantity"]


def window_text(window):
    return f"{window['quantity']} {window['low']:g} to {window['high']:g} (falls off over {window['sigma_low']:g} below, {window['sigma_high']:g} above)"


def energy_text(energy):
    if energy["source"] == "light":
        return f"light, {energy['band_nm'][0]:g} to {energy['band_nm'][1]:g} nm, efficiency {energy['efficiency']:g}"
    return "chemical (internal and tidal heat)"


def supplies_of(row):
    """The supplies of a class as a list (the rule file holds one object or a list of them)."""
    supply = row.get("supply")
    return [] if not supply else supply if isinstance(supply, list) else [supply]


def supply_text(supplies):
    if not supplies:
        return "none"
    return "; ".join(f"{supply['quantity']} (full at {supply['ref']:g} of Earth's, placeholder {supply['placeholder']:g})" for supply in supplies)


def life_text(life):
    oxidant_yield = "from the air" if life["oxidant_yield"] is None else f"{life['oxidant_yield']:g}"
    return f"{life['metabolism']}; breathes {life['oxidant']} (yield {oxidant_yield}); first step x{life['first_step_scale']:g}"


def render(document):
    lines = [
        "# Life class rules", "",
        "Generated from `src/planetgen/physics/data/life_classes.json` by `scripts/build_life_class_rules.py`; do not edit by hand. "
        "The classifier (`planetgen.physics.life_classes`, GEN.886) reads that file and nothing else; the design is "
        "[life-classes.md](life-classes.md). Every number in the file carries a source tier (T1 to T9) and a tweak level (L0 to L4) "
        "(see [stand-in-values.md](../research/stand-in-values.md)).", "",
        f"Rule file version {document['version']}. Reference star {document['reference_star_k']:g} K, reference band "
        f"{document['reference_band_nm'][0]:g} to {document['reference_band_nm'][1]:g} nm.", "",
        "## The classes", "",
        "| Class | Name | Drinks | Eats | Breathes | Gives off | Reaches | Tier |", "|---|---|---|---|---|---|---|---|",
    ]
    for row in document["classes"]:
        lines.append("| " + " | ".join(_cell(row[key]) for key in ("id", "name", "drinks", "eats", "breathes", "gives_off", "reaches", "tier")) + " |")
    lines += ["", "## What each class needs", "",
              "| Class | Solvent | Energy | Gates (all must hold) | Windows | Supply |", "|---|---|---|---|---|---|"]
    for row in document["classes"]:
        gates = "; ".join(gate_text(g) for g in row["gates"]) or "none"
        windows = "; ".join(window_text(w) for w in row["windows"]) or "none"
        solvent = f"{row['solvent']['medium']} (factor {row['solvent']['factor']:g})"
        lines.append(f"| {row['id']} | {_cell(solvent)} | {_cell(energy_text(row['energy']))} | {_cell(gates)} | {_cell(windows)} | "
                     f"{_cell(supply_text(supplies_of(row)))} |")
    lines += ["", "## What each class lives on", "", "| Class | Metabolism, oxidant and first step | Colour in visible light | Outside visible light |", "|---|---|---|---|"]
    for row in document["classes"]:
        life = row["life"]
        lines.append(f"| {row['id']} | {_cell(life_text(life))} | {_cell(', '.join(life['reflection_spectrum_visible']) or 'none')} | "
                     f"{_cell(', '.join(life['reflection_spectrum_non_visible']) or 'none')} |")
    lines += ["", "## Why each number is what it is", ""]
    for row in document["classes"]:
        lines += [f"### Class {row['id']}: {row['name']}", ""]
        for label, entry in [("solvent", row["solvent"]), ("energy", row["energy"])] + [("gate", g) for g in row["gates"]] \
                + [("window", w) for w in row["windows"]] + [("supply", s) for s in supplies_of(row)] + [("life", row["life"])]:
            what = {"solvent": row["solvent"]["medium"], "energy": row["energy"]["source"], "life": row["life"]["metabolism"]}.get(label) or entry.get("quantity", "")
            lines.append(f"- {label} {what}: {entry['tier']}, {entry['tweak']}. {entry['basis']}")
        if row.get("note"):
            lines.append(f"- note: {row['note']}")
        lines.append("")
    chem = document["non_water_chemistry"]
    lines += ["## Classes that do not drink water", "",
              f"The hydrogen inventory of the chemistry term is held at {chem['hydrogen_inventory_floor']:g} Earth units: {chem['tier']}, {chem['tweak']}. {chem['basis']}.", ""]
    return "\n".join(lines)


def load():
    with open(RULES, encoding="utf-8") as handle:
        return json.load(handle)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    text = render(load())
    if args.check:
        with open(OUT, encoding="utf-8") as handle:
            if handle.read() != text:
                print("docs/design/habitability/life-class-rules.md is out of date: run python scripts/build_life_class_rules.py")
                return 1
        return 0
    with open(OUT, "w", encoding="utf-8") as handle:
        handle.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
