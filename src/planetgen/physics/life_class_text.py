# planetgen/physics/life_class_text.py

"""
The words of the life classes (GEN.892, docs/design/habitability/life-classes.md sections 4 and 5)
===============================================================================================

What a class drinks, eats, breathes and gives off, and what each rung of the stage ladder is called for it, in plain words for the
planet page, the search results and the API. The lines of the first kind are the rule file's own (`physics/data/life_classes.json`);
the stage words are `physics/data/life_class_text.json`:

* `rungs`: a name and a sentence for every rung of the ladder, for any class (`{n}` is the class's whole id when one is given, "3a", else its number);
* `grades`: the three animal grades A1 to A3 of rung 5;
* `classes`: a name and a sentence for the rungs where a class has words of its own. A class has words for the rungs the design says
  it may reach; any other rung falls back to the generic words, so a strange-solvent class that does reach complex cells in some
  draw still reads correctly. Nothing here says a rung is out of reach: the physics decides that (`life_classes.py`).

Use:

    stage("1", 2)                    # ("Oxygenating mats", "Light-eating microbes spread ...")
    stage("4a", 5, grade="A2")       # the words of the animal grade, with the class's own words for rung 5 kept as its sentence
    lines("7")                       # {"drinks": ..., "eats": ..., "breathes": ..., "gives_off": ...}
    sentence("1a")                   # "It drinks water, eats light 400 to 700 nm ..., breathes ... and gives off ..."
"""

import json
import os
import re

from planetgen.physics import elements, life_classes

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "life_class_text.json")

RUNGS = (0, 1, 2, 3, 4, 5, 6)
GRADES = ("A1", "A2", "A3")


class LifeClassTextError(ValueError):
    """The text file is malformed, or a class, rung or grade asked for does not exist."""


def _check(document):
    for key in ("rungs", "grades", "classes"):
        if key not in document:
            raise LifeClassTextError(f"life_class_text.json has no {key!r}")
    if sorted(document["rungs"]) != [str(r) for r in RUNGS]:
        raise LifeClassTextError("life_class_text.json needs a generic entry for each rung 0 to 6")
    if sorted(document["grades"]) != sorted(GRADES):
        raise LifeClassTextError("life_class_text.json needs an entry for each animal grade A1 to A3")
    numbers = {str(c.number) for c in life_classes.CLASSES}
    for number, rungs in document["classes"].items():
        if number not in numbers:
            raise LifeClassTextError(f"words for class {number}, which is not a life class")
        for rung in rungs:
            if rung not in {str(r) for r in RUNGS}:
                raise LifeClassTextError(f"class {number} has words for rung {rung}, which is not a rung")
    for entry in [*document["rungs"].values(), *document["grades"].values(),
                  *(e for rungs in document["classes"].values() for e in rungs.values())]:
        if not (entry.get("name") and entry.get("description")):
            raise LifeClassTextError(f"an entry needs a name and a description: {entry!r}")
    missing = numbers - set(document["classes"])
    if missing:
        raise LifeClassTextError(f"no words for class {sorted(missing)}")


def _load():
    with open(DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    _check(document)
    return document


DOCUMENT = _load()
VERSION = DOCUMENT["version"]


def number_of(class_id):
    """The class number of `class_id` (`"4c"`, `4`, `"4"` all give 4)."""
    if isinstance(class_id, int):
        number = class_id
    else:
        match = re.fullmatch(r"(\d+)[a-z]?", str(class_id))
        if not match:
            raise LifeClassTextError(f"{class_id!r} is not a life class")
        number = int(match.group(1))
    if str(number) not in DOCUMENT["classes"]:
        raise LifeClassTextError(f"{class_id!r} is not a life class")
    return number


def stage(class_id, rung, grade=None):
    """
    The words of one rung of the stage ladder for a class.

    Parameters:
        class_id: A class number or id (`"1a"`, `4`).
        rung: 0 to 6.
        grade: For rung 5, `"A1"`, `"A2"` or `"A3"`: the name is then the grade's and the description is the rung's sentence for the class followed by the grade's.

    Returns:
        tuple: `(name, description)`.
    """
    number = number_of(class_id)
    name = str(class_id) if str(class_id) in life_classes.BY_ID else str(number)
    if rung not in RUNGS:
        raise LifeClassTextError(f"{rung!r} is not a rung (0 to 6)")
    if grade is not None:
        if rung != 5:
            raise LifeClassTextError("only rung 5 (animals) has grades")
        if grade not in GRADES:
            raise LifeClassTextError(f"{grade!r} is not an animal grade (A1 to A3)")
        entry = DOCUMENT["grades"][grade]
        own = DOCUMENT["classes"][str(number)].get("5") or DOCUMENT["rungs"]["5"]
        return entry["name"], own["description"].replace("{n}", name) + " " + entry["description"]
    entry = DOCUMENT["classes"][str(number)].get(str(rung)) or DOCUMENT["rungs"][str(rung)]
    return entry["name"], entry["description"].replace("{n}", name)


def stage_name(class_id, rung, grade=None):
    return stage(class_id, rung, grade)[0]


def stage_description(class_id, rung, grade=None):
    return stage(class_id, rung, grade)[1]


def lines(class_id):
    """What the class drinks, eats, breathes and gives off, from the rule file: `{"drinks": ..., "eats": ..., "breathes": ..., "gives_off": ...}`.
    A bare class number (`"4"`) reads the first variant of it (`"4a"`)."""
    row = life_classes.BY_ID.get(str(class_id)) or life_classes.describe(_first_variant(class_id))
    return {key: getattr(row, key) for key in ("drinks", "eats", "breathes", "gives_off")}


def _first_variant(class_id):
    number = number_of(class_id)
    for row in life_classes.CLASSES:
        if row.number == number:
            return row.id
    raise LifeClassTextError(f"{class_id!r} is not a life class")


def sentence(class_id):
    """The four lines as one plain sentence: "It drinks ..., eats ..., breathes ... and gives off ...." """
    said = lines(class_id)
    return f"It drinks {said['drinks']}, eats {said['eats']}, breathes {said['breathes']} and gives off {said['gives_off']}."
