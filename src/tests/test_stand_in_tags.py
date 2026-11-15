"""
The stand-in tags of the code and the registry they name (GEN.1632).

A tag in the code (`[Stand-in 6.294 AD1, tier T7, tweak L1.]`) names the section of the registry beside the row id, because the same id is used in
several sections of `docs/design/research/gap-filling-plan.md` (an id is unique in its section, not in the plan). The section and the id name one
row; the tier and the tweak the tag states are the row's. Every table row of the plan that carries a tier and a tweak is in the registry, and a
row that names a `module.CONSTANT` names one that exists.
"""

import importlib
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts"))
import build_stand_in_values as registry  # noqa: E402

SOURCE = os.path.join(HERE, "..", "planetgen")
TAG = re.compile(r"\[Stand-in (?P<section>\d+\.\d+[a-z]?) (?P<row>[A-Z]+\d+[a-z]?)(?P<rest>[^\]]*)\]")
ANY_TAG = re.compile(r"\[Stand-in\b")
STATED = re.compile(r"^, tier (?P<tier>T\d)(?:[^,;]*?)?, tweak (?:level )?(?P<tweak>L\d)")
CONSTANT = re.compile(r"\b([a-z][a-z0-9_]*)\.([A-Z][A-Z0-9_]{1,}[A-Z0-9])\b")


def _plan():
    with open(registry.PLAN, encoding="utf-8") as stream:
        return stream.read()


def _rows():
    """`{(section number, row id): (tier code, tweak, constant cell)}` of the registry, and the ids a section repeats."""
    tables = registry.registry_tables(_plan())
    numbers = {key: number for key, number, _, _ in registry.sections()}
    rows, repeated = {}, []
    for key, (header, body) in tables.items():
        lowered = [name.lower() for name in header]
        for row in body:
            at = (numbers[key], row[0])
            if at in rows:
                repeated.append(f"{at[0]} {at[1]}")
            constant = row[lowered.index("constant")] if "constant" in lowered else row[lowered.index("stand-in")] if "stand-in" in lowered else " ".join(row[4:5])
            rows[at] = (re.match(r"T\d", row[lowered.index("tier")]).group(0), row[lowered.index("tweak")], constant)
    return rows, repeated


def _tags():
    """`[(file, line, text)]` of every `[Stand-in ...]` tag in the code, one per tag."""
    found = []
    for root, _, names in os.walk(SOURCE):
        for name in names:
            if not name.endswith(".py"):
                continue
            path = os.path.join(root, name)
            with open(path, encoding="utf-8") as stream:
                text = stream.read()
            for match in ANY_TAG.finditer(text):
                end = text.index("]", match.start())
                found.append((os.path.relpath(path, SOURCE), text.count("\n", 0, match.start()) + 1, text[match.start():end + 1]))
    return found


def test_a_row_id_is_not_repeated_inside_a_section():
    """The section and the id name one row: no section holds the same id twice."""
    assert _rows()[1] == []


def test_every_registry_table_row_of_the_plan_is_in_the_registry():
    """A row of a table the registry does not read is a stand-in nothing lists: no table is left out, none ends at a blank line."""
    tables = registry.registry_tables(_plan())
    captured = {tuple(row) for _, body in tables.values() for row in body}
    row = re.compile(r"\| [A-Za-z0-9`:.\-]+ \| [A-Z] \| T\d[^|]*\| L\d[^|]*\|")
    left = [number for number, line in enumerate(_plan().splitlines(), 1) if row.match(line) and tuple(registry._cells(line)) not in captured]
    assert left == [], f"rows of the plan outside the registry tables, by line: {left}"


def test_every_tag_of_the_code_names_its_section_and_an_id():
    """`[Stand-in 6.294 AD1, ...]`: the section number first, since an id alone is used by several sections."""
    bad = [f"{name}:{line}: {text[:60]!r}" for name, line, text in _tags() if not TAG.fullmatch(text.replace("\n", " "))]
    assert bad == [], "tags with no section and id:\n" + "\n".join(bad)


def test_every_tag_of_the_code_names_a_row_of_the_registry():
    rows = _rows()[0]
    bad = [f"{name}:{line}: {match['section']} {match['row']}" for name, line, text in _tags()
           if (match := TAG.fullmatch(text.replace("\n", " "))) and (match["section"], match["row"]) not in rows]
    assert bad == [], "tags that name no registry row:\n" + "\n".join(bad)


def test_the_tier_and_the_tweak_a_tag_states_are_the_rows():
    """Where the tag says `tier T6, tweak L1` straight after its id, the row has that tier and that tweak."""
    rows = _rows()[0]
    bad = []
    for name, line, text in _tags():
        match = TAG.fullmatch(text.replace("\n", " "))
        if not match or (match["section"], match["row"]) not in rows:
            continue
        stated = STATED.match(match["rest"])
        if stated is None:
            continue
        tier, tweak, _ = rows[(match["section"], match["row"])]
        if (stated["tier"], stated["tweak"]) != (tier, tweak):
            bad.append(f"{name}:{line}: {match['section']} {match['row']} says {stated['tier']} {stated['tweak']}, the registry {tier} {tweak}")
    assert bad == [], "tags that differ from their row:\n" + "\n".join(bad)


def test_a_row_names_only_constants_that_exist():
    """`module.CONSTANT` in the constant cell of a row is a name the code has: a row of a constant that was removed is retired."""
    missing = []
    for (section, row), (_, _, cell) in sorted(_rows()[0].items()):
        for name, constant in CONSTANT.findall(cell):
            found = _module(name)
            if found is not None and not hasattr(found, constant):
                missing.append(f"{section} {row}: {name}.{constant}")
    assert missing == [], "rows that name constants that no longer exist:\n" + "\n".join(missing)


def _module(name):
    """The module `planetgen.<package>.<name>` where one file of that name is unique in the code, else `None`."""
    paths = [root for root, _, names in os.walk(SOURCE) if f"{name}.py" in names]
    if len(paths) != 1:
        return None
    package = os.path.relpath(paths[0], os.path.join(SOURCE, "..")).replace(os.sep, ".")
    try:
        return importlib.import_module(f"{package}.{name}")
    except Exception:
        return None
