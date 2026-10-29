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


def supply_text(supply):
    if not supply:
        return "none"
    return f"{supply['quantity']} (placeholder {supply['placeholder']:g} until it is stored)"


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
                     f"{_cell(supply_text(row.get('supply')))} |")
    lines += ["", "## Why each number is what it is", ""]
    for row in document["classes"]:
        lines += [f"### Class {row['id']}: {row['name']}", ""]
        for label, entry in [("solvent", row["solvent"]), ("energy", row["energy"])] + [("gate", g) for g in row["gates"]] \
                + [("window", w) for w in row["windows"]] + ([("supply", row["supply"])] if row.get("supply") else []):
            what = {"solvent": row["solvent"]["medium"], "energy": row["energy"]["source"]}.get(label) or entry.get("quantity", "")
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
