#!/usr/bin/env python3
# scripts/build_atmosphere_recipes.py

"""
Builds src/planetgen/physics/data/atmosphere_recipes.json (GEN.478)
===================================================================

The 104 rows of part 09 of the atmosphere research (docs/design/composition/atmosphere-composition-parts/09-weighted-composition-tables.md)
in 16 tables a to p. The words of each row (option, gate, weight, weight source, status, notes) are read from the part itself;
the machine columns (kind, classes, conditions, weight expression, distribution, constants, corrections) are in
`atmosphere_recipes_data.py`. The result is written for the loader `planetgen.physics.atmosphere_recipes`.

    python3 scripts/build_atmosphere_recipes.py            write the file
    python3 scripts/build_atmosphere_recipes.py --check    fail if the file is not what this would write

A row's source state follows its weight source: a row whose source cell carries only READ, COMPUTED or DATA labels is
`verified`; one that also carries 2ND-HAND, SNIPPET, RECALLED, NOT REACHED, UNSUPPORTED or DESIGN labels is `partly verified`
with those labels as its limit (a derived row whose edge width is a design choice is usable, and says so); a row with no
sourced weight (status DESIGN) is `open` and names the research item (R1 to R15 = GEN.480 to GEN.494) that replaces it; a row
built from a measured Solar System body is a `stand-in` with its replacement item (Boss, 1:44 AM EDT 2026-10-17).
"""

import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import atmosphere_recipes_data as data  # noqa: E402

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
PART = os.path.join("docs", "design", "composition", "atmosphere-composition-parts", "09-weighted-composition-tables.md")
OUT = os.path.join(ROOT, "src", "planetgen", "physics", "data", "atmosphere_recipes.json")

CATEGORIES = {}
LABELS = ("READ", "2ND-HAND", "DATA", "COMPUTED", "SNIPPET", "RECALLED", "NOT REACHED", "UNSUPPORTED", "DESIGN")
STRONG = {"READ", "COMPUTED", "DATA"}
STATUSES = ("SOURCED", "DERIVED", "DESIGN")
RETRIEVED = "2026-10-18"


def parse_part():
    """The rows of the part's section 2: [(category, cells)], and the category titles."""
    with open(os.path.join(ROOT, PART), encoding="utf-8") as handle:
        text = handle.read()
    rows, category = [], None
    for line in text.splitlines():
        heading = re.match(r"### ([a-p])\. (.*)", line)
        if heading:
            category = heading.group(1)
            CATEGORIES[category] = heading.group(2).strip()
            continue
        if line.startswith("## 3."):
            category = None
        if category and line.startswith("| ") and not line.startswith("| id") and not line.startswith("|---"):
            cells = [c.strip() for c in line.strip().strip("|").split(" | ")]
            if len(cells) != 7:
                raise SystemExit(f"row of {category} has {len(cells)} cells: {line[:80]}")
            rows.append((category, cells))
    return rows


def labels_of(source):
    found = []
    for label in LABELS:
        if re.search(r"(?<![A-Z0-9-])" + re.escape(label) + r"(?![A-Z0-9-])", source):
            found.append(label)
    return found


def state_of(entry):
    """(state, limit, replacement, research items) of a row."""
    status, labels = entry["status"], entry["labels"]
    if entry.get("state_override") == "stand-in":
        return "stand-in", "", entry["replacement"], entry["research"]
    if entry.get("state_override") == "partly verified":
        return "partly verified", entry["override_limit"], "", []
    if status == "DESIGN":
        return "open", "", "", entry["research"]
    weak = [label for label in labels if label not in STRONG]
    if weak:
        return "partly verified", "part of the weight source is " + ", ".join(weak) + " (its edge widths or constants are design choices)", "", []
    return "verified", "", "", []


def build():
    rows = parse_part()
    ids = [cells[0] for _cat, cells in rows]
    if len(ids) != 104 or len(set(ids)) != 104:
        raise SystemExit(f"part 09 has {len(ids)} rows ({len(set(ids))} distinct), not 104")
    if set(ids) != set(data.MACHINE):
        raise SystemExit("the machine columns and the part disagree: "
                         f"missing {sorted(set(ids) - set(data.MACHINE))}, extra {sorted(set(data.MACHINE) - set(ids))}")
    counts = {s: sum(1 for _c, cells in rows if cells[5].split()[0] == s) for s in STATUSES}
    if counts != {"SOURCED": 8, "DERIVED": 67, "DESIGN": 29}:
        raise SystemExit(f"part 09 counts changed: {counts}")
    rows = [(category, cells) for category, cells in rows if cells[0] not in data.RETIRED]      # GEN.524 part 2: the part keeps them as the record
    entries = []
    resolved_source = {}
    previous = None
    for category, cells in rows:
        row_id, option, gate, weight, source, status, notes = cells
        machine = dict(data.MACHINE[row_id])
        parts = status.split()
        research = [r for r in parts[1:] if re.fullmatch(r"R\d+", r)]
        status_word = machine.pop("state", None)
        if status_word in ("derived",):
            status, research = "DERIVED", []
        entry = {
            "id": row_id,
            "category": category,
            "kind": machine["kind"],
            "option": option,
            "gate_text": gate,
            "weight_text": weight,
            "source": machine.pop("source", source),
            "status": status.split()[0],
            "research": research,
            "notes": notes,
        }
        if status_word == "stand-in":
            entry["state_override"] = "stand-in"
            entry["replacement"] = machine.pop("replacement")
        elif status_word == "partly verified":
            entry["state_override"] = "partly verified"
            entry["override_limit"] = machine.pop("limit")
        correction = machine.pop("correction", "")
        if correction:
            entry["notes"] = (notes + " " + correction).strip()
        text = entry["source"]
        if re.match(r"same\b", text) and previous is not None:
            text = resolved_source[previous]
        resolved_source[row_id] = text
        previous = row_id
        entry["labels"] = labels_of(text)
        state, limit, replacement, items = state_of(entry)
        entry.pop("state_override", None)
        entry.pop("override_limit", None)
        entry["state"] = state
        entry["limit"] = limit
        entry["replacement"] = replacement
        entry["items"] = [data.RESEARCH[r] for r in research] if state == "open" else (
            [replacement] if replacement else [])
        machine.pop("replacement", None)
        entry["classes"] = machine.pop("classes")
        entry["conditions"] = machine.pop("when", [])
        for key in list(machine):
            entry[key] = machine.pop(key)
        entries.append(entry)
    groups = {}
    for name, spec in data.GROUPS.items():
        groups[name] = {"rows": spec["rows"], "fallback": spec["fallback"], "constants": spec.get("constants", {})}
    for name, spec in groups.items():
        for row_id in spec["rows"]:
            owner = next(e for e in entries if e["id"] == row_id)
            if owner.get("kind") != "option" or owner.get("group") != name:
                raise SystemExit(f"row {row_id} is listed in group {name} but is {owner.get('kind')} of {owner.get('group')}")
        if spec["fallback"] not in (None, "none") and spec["fallback"] not in spec["rows"]:
            raise SystemExit(f"group {name}: the fallback {spec['fallback']} is not one of its rows")
    for entry in entries:
        if entry["kind"] == "option" and entry["id"] not in groups[entry["group"]]["rows"]:
            raise SystemExit(f"option {entry['id']} is in no group")
    return {
        "title": "The atmosphere possibility tables (GEN.478): 16 categories, 104 rows",
        "retrieved": RETRIEVED,
        "part": PART.replace(os.sep, "/"),
        "version": 1,
        "categories": CATEGORIES,
        "constants": data.CONSTANTS,
        "research": data.RESEARCH,
        "groups": groups,
        "entries": entries,
    }


def render(document):
    return json.dumps(document, indent=1, ensure_ascii=False) + "\n"


def main(argv):
    text = render(build())
    if "--check" in argv:
        with open(OUT, encoding="utf-8") as handle:
            if handle.read() != text:
                print(f"{OUT} is not what this builder writes; run it without --check", file=sys.stderr)
                return 1
        print("atmosphere_recipes.json matches the builder")
        return 0
    with open(OUT, "w", encoding="utf-8") as handle:
        handle.write(text)
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
