"""
Typed commands (ADM.121): a line such as `add sectors around sector 9 radius 20` becomes the Add sectors form filled in
and opened at its Review step. Nothing runs from a command: the Review step is the same one the form leads to, and Run is
still a button on it. The grammar is small on purpose.

    add sectors [radius N[pc|ly]] [density N]                    around a random populated spot
    add sectors around sector N radius R                         around a filled sector
    add sectors around shell S [band B] wedge W radius R         around a sector address
    add sectors around X Y [Z] radius R                          around a position in parsecs
    add sectors shell S [band B] [limit N]                       a whole shell
    add sectors shells A:B [bands A:B] [wedges A:B] [limit N]    a range of the grid
    check [deep]                                                  check the database
    scatter | new galaxy | wipe                                  open that task's form

`parse` returns `None` for a line that is not a command, a `Command` for one that is, and raises `CommandError` (with a
message for the admin) for a line that starts like a command but cannot be read.
"""

import math
from dataclasses import dataclass, field

from planetgen.physics.units import ly_to_pc

VERBS = (
    ("add sectors", "add-sectors"), ("add", "add-sectors"), ("sectors", "add-sectors"),
    ("scatters", "scatters"), ("scatter", "scatters"),
    ("new galaxy", "new-galaxy"), ("plan only", "plan-only"), ("plan", "plan-only"), ("wipe", "wipe"), ("check", "check"),
)
"""tuple: The words that start a command, longest first, and the task each opens."""

OPEN_ONLY = frozenset({"scatters", "new-galaxy", "plan-only", "wipe"})
"""frozenset: Tasks a command only opens the form of: they need choices or a typed confirmation a line cannot carry."""

VALUE_WORDS = {
    "radius": "radius", "density": "density", "shell": "shell", "band": "band", "wedge": "wedge", "sector": "sector",
    "shells": "shells", "bands": "bands", "wedges": "wedges", "limit": "limit",
}


FILLER = frozenset({"around", "at", "near", "of", "within"})
"""frozenset: Words that only make a line read well; they are ignored."""


class CommandError(ValueError):
    """A line that starts like a command but cannot be read; the message says why."""


@dataclass(frozen=True)
class Command:
    """A parsed command: the task it opens, the form fields it fills, whether to open the Review step, and a sentence
    saying what it understood."""

    task: str
    fields: dict = field(default_factory=dict)
    review: bool = True
    summary: str = ""


def _number(word, label, kind=float):
    try:
        value = kind(word)
    except ValueError:
        raise CommandError(f"{label} needs a number, not '{word}'.") from None
    if isinstance(value, float) and not math.isfinite(value):
        raise CommandError(f"{label} needs a finite number.")
    return value


def _radius_pc(word):
    """A radius word: parsecs, or light-years with an `ly` suffix (`40ly`), in parsecs rounded to 0.1."""
    text = word.lower()
    if text.endswith("ly"):
        return str(round(ly_to_pc(_number(text[:-2], "Radius")), 1))
    return str(_number(text[:-2] if text.endswith("pc") else text, "Radius"))


def _split(text):
    """`(task, rest words)` for a line that starts with a verb, else `(None, [])`."""
    words = text.lower().split()
    line = " ".join(words)
    for verb, task in VERBS:
        if line == verb or line.startswith(verb + " "):
            return task, line[len(verb):].split()
    return None, []


def parse(text):
    """The `Command` for a typed line, `None` when it is not one, or `CommandError`."""
    task, words = _split(text)
    if task is None:
        return None
    if task == "check":
        deep = [w for w in words if w != "deep"]
        if deep:
            raise CommandError(f"Check takes only 'deep', not '{deep[0]}'.")
        return Command("check", {"action": "check_db", **({"deep": "1"} if words else {})}, True,
                       "Check the database" + (", deeply" if words else ""))
    if task in OPEN_ONLY:
        if words:
            raise CommandError("That task opens its form; choose its options there.")
        return Command(task, {}, False, "Open the form")
    return _add_sectors(words)


def _add_sectors(words):
    values, position = {}, []
    words = [w for w in words if w not in FILLER]
    index = 0
    while index < len(words):
        word = words[index]
        key = VALUE_WORDS.get(word)
        if key is None:
            if _is_number(word):
                position.append(word)
                index += 1
                continue
            raise CommandError(f"I don't know '{word}'. Try 'radius 12', 'sector 9', 'shell 40' or 'shells 10:20'.")
        if index + 1 >= len(words):
            raise CommandError(f"'{word}' needs a value after it.")
        if key in values:
            raise CommandError(f"'{word}' is given twice.")
        values[key] = words[index + 1]
        index += 2
    return _build(values, position)


def _is_number(word):
    try:
        float(word)
    except ValueError:
        return False
    return True


def _build(values, position):
    grid = {"shells", "bands", "wedges"} & values.keys()
    if position and len(position) not in (2, 3):
        raise CommandError("A position is x y or x y z, in parsecs.")
    if position and (values.keys() - {"radius"}):
        raise CommandError("A position cannot be combined with a shell, band, wedge or sector.")
    if "sector" in values and (values.keys() & {"shell", "band", "wedge"} or grid):
        raise CommandError("A sector cannot be combined with a shell, band or wedge.")
    if grid and (values.keys() & {"shell", "band", "wedge", "sector", "radius", "density"}):
        raise CommandError("A range (shells, bands, wedges) cannot be combined with a place or a radius.")
    if position:
        x, y, *z = position
        fields = {"what": "place", "place": "position", "center_x_pc": x, "center_y_pc": y,
                  "center_radius_pc": _need_radius(values)}
        if z:
            fields["center_z_pc"] = z[0]
        return Command("add-sectors", fields, True, f"Add sectors within {fields['center_radius_pc']} pc of {', '.join(position)} pc")
    if "sector" in values:
        _number(values["sector"], "Sector", int)
        fields = {"what": "place", "place": "sector", "center_sector": values["sector"],
                  "center_radius_pc": _need_radius(values)}
        return Command("add-sectors", fields, True,
                       f"Add sectors within {fields['center_radius_pc']} pc of sector {values['sector']}")
    if grid:
        fields = {"what": "grid", "mode": "span"}
        for name, field_name in (("shells", "span_shells"), ("bands", "span_bands"), ("wedges", "span_wedges")):
            if name in values:
                fields[field_name] = values[name]
        if "limit" in values:
            _number(values["limit"], "Limit", int)
            fields["span_limit"] = values["limit"]
        return Command("add-sectors", fields, True, "Add sectors across a range of the grid")
    if "wedge" in values:
        if "shell" not in values:
            raise CommandError("A wedge needs its shell: 'shell 40 wedge 7 radius 12'.")
        fields = {"what": "place", "place": "address", "center_shell": values["shell"], "center_wedge": values["wedge"],
                  "center_radius_pc": _need_radius(values)}
        if "band" in values:
            fields["center_band"] = values["band"]
        return Command("add-sectors", fields, True,
                       f"Add sectors within {fields['center_radius_pc']} pc of shell {values['shell']} wedge {values['wedge']}")
    if "shell" in values:
        _number(values["shell"], "Shell", int)
        fields = {"what": "grid", "mode": "shell", "shell_index": values["shell"]}
        if "band" in values:
            fields["shell_band"] = values["band"]
        if "limit" in values:
            _number(values["limit"], "Limit", int)
            fields["shell_limit"] = values["limit"]
        return Command("add-sectors", fields, True, f"Add the sectors of shell {values['shell']}")
    if values.keys() & {"band", "limit"}:
        raise CommandError("A band or limit goes with a shell: 'shell 40 band 1 limit 200'.")
    fields = {"what": "place", "place": "random"}
    if "radius" in values:
        fields["radius_pc"] = _radius_pc(values["radius"])
    if "density" in values:
        fields["min_start_density"] = str(_number(values["density"], "Density"))
    return Command("add-sectors", fields, True, "Add sectors around a random populated spot")


def _need_radius(values):
    if "radius" not in values:
        raise CommandError("This place needs a radius: add 'radius 12'.")
    return _radius_pc(values["radius"])
