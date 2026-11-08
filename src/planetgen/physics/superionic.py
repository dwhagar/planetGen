# planetgen/physics/superionic.py

"""
GEN.1422: where water and ammonia are superionic. A deep layer of an ice giant or a water world holds water and ammonia at the pressure and temperature of its
bottom; this module says whether they are there

* an ice (or a molecular fluid) under the floor,
* **superionic** between the floor and the ceiling: a solid of oxygen (or nitrogen) atoms whose hydrogen ions move freely (bcc, then fcc superionic ice,
  [R5171]; the phase of Millot's ice XVIII, [R3081]; the superionic ammonia of [R5172]),
* an **ionic fluid** above the ceiling (the melt line), and
* a **metallic fluid** above the temperature where the fluid turns metallic (7000 K for water, 5500 K for ammonia, [R5172]).

The floor and the ceiling of each species are tables of `(pressure GPa, temperature K)` points in `physics/data/superionic_boundary.json`, with their sources. Between
points the line is straight; past the last point the last temperature holds; under the first point there is no superionic phase of that species. Where the
experiments disagree on where the phase starts the lowest floor is used (Boss's rule of the widest range). It replaces the single pressure of 100 GPa that
`giant_interior` held.

Use:

    superionic.state("water", 80.0, 3000.0)      # "superionic"
    superionic.is_superionic("ammonia", 20.0, 2000.0)
"""

import json
import math
import os

from planetgen.physics import elements
from planetgen.physics.elements import DataError

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "superionic_boundary.json")
"""str: The data file."""

SPECIES = ("water", "ammonia")
"""tuple: The species the table covers."""

STATES = ("solid_or_molecular", "superionic", "ionic_fluid", "metallic_fluid")
"""tuple: What `state` answers."""

FALLBACK_FLOOR_BEYOND_GPA = 60.0
"""float: STAND-IN SI1 (T6, L0): the water floor measured up to this pressure (1900 K, [R5171]) is held at all higher pressures (the table's last floor point)."""

CEILING_INTERPOLATION_GPA = (20.0, 105.0)
"""tuple: STAND-IN SI2 (T6, L1): the water ceiling is the straight line between its table points at these pressures (the triple point at 20 GPa and the point at 105 GPa, [R5171])."""

AMMONIA_FLOOR_GPA_K = (30.0, 1500.0)
"""tuple: STAND-IN SI3 (T7, L1): ammonia is superionic from this pressure and temperature upward (recalled, not opened: [R5172] gives the range only)."""

FALLBACK_TEMPERATURE_K = 2500.0
"""float: STAND-IN SI4 (T7, L1): the temperature a layer is taken to have when its structure model carries none (the constant-density
stand-in used when the hydrogen-helium tables are missing): a typical temperature of an ice giant's deep ice layer on the planetary isentrope of [R5172], 2000 to 5000 K."""

_cache = {}


def _curve(raw, what):
    points = []
    for item in raw or []:
        if len(item) != 2 or not all(isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) and v > 0 for v in item):
            raise DataError(f"{what} must be pairs of positive numbers, not {item!r}")
        points.append((float(item[0]), float(item[1])))
    if not points or any(b[0] <= a[0] for a, b in zip(points, points[1:])):
        raise DataError(f"{what} needs points with rising pressure")
    return tuple(points)


def load(path=DATA_FILE):
    """Reads and checks the table: `{species: {"material", "floor", "ceiling", "metallic_k"}}`; a malformed row is refused (`DataError`)."""
    try:
        with open(path, encoding="utf-8") as handle:
            document = json.load(handle)
    except (OSError, ValueError) as error:
        raise DataError(f"cannot read the superionic boundary table {path}: {error}") from error
    table = {}
    for name in SPECIES:
        raw = (document.get("species") or {}).get(name)
        if not raw:
            raise DataError(f"the superionic boundary table has no {name} row")
        metallic = raw.get("metallic_k")
        if not isinstance(metallic, (int, float)) or metallic <= 0:
            raise DataError(f"{name} needs a metallic temperature")
        table[name] = {"material": raw.get("material", ""), "floor": _curve(raw.get("floor"), f"the {name} floor"),
                       "ceiling": _curve(raw.get("ceiling"), f"the {name} ceiling"), "metallic_k": float(metallic)}
    water, ammonia = table["water"], table["ammonia"]
    if (water["floor"][-1][0] != FALLBACK_FLOOR_BEYOND_GPA or (water["ceiling"][0][0], water["ceiling"][1][0]) != CEILING_INTERPOLATION_GPA
            or ammonia["floor"][0] != AMMONIA_FLOOR_GPA_K):
        raise DataError("the stand-in points SI1 to SI3 (section 6.233) do not match the table's rows")
    return table


def table():
    """The checked table, read once."""
    if "table" not in _cache:
        _cache["table"] = load()
    return _cache["table"]


def _along(points, pressure_gpa):
    """The temperature of the line through `points` at `pressure_gpa`: `None` under the first point, the last temperature past the last."""
    if pressure_gpa < points[0][0]:
        return None
    for (p0, t0), (p1, t1) in zip(points, points[1:]):
        if pressure_gpa <= p1:
            return t0 + (t1 - t0) * (pressure_gpa - p0) / (p1 - p0)
    return points[-1][1]


def floor_k(species, pressure_gpa):
    """The temperature (K) where `species` turns superionic at `pressure_gpa`; `None` where it is not superionic at any temperature (under its first point)."""
    return _along(table()[species]["floor"], pressure_gpa)


def ceiling_k(species, pressure_gpa):
    """The temperature (K) where superionic `species` melts to an ionic fluid at `pressure_gpa` (`None` under the first point)."""
    return _along(table()[species]["ceiling"], pressure_gpa)


def state(species, pressure_gpa, temperature_k):
    """What `species` is at the pressure and temperature: one of `STATES`."""
    row = table()[species]
    if temperature_k >= row["metallic_k"]:
        return "metallic_fluid"
    floor, ceiling = floor_k(species, pressure_gpa), ceiling_k(species, pressure_gpa)
    if floor is None or temperature_k < floor:
        return "solid_or_molecular"
    return "superionic" if temperature_k <= ceiling else "ionic_fluid"


def is_superionic(species, pressure_gpa, temperature_k):
    """Whether `species` is a superionic solid there."""
    return state(species, pressure_gpa, temperature_k) == "superionic"


def layer_state(components, pressure_gpa, temperature_k):
    """The state of a layer of water and ammonia (the `components` mass mapping of a `giant_interior.Layer`) at its bottom: the state of the species that holds the most
    mass of the two; `None` when the layer holds neither. A layer of both is given the state of the larger share."""
    water = (components.get("water_liquid", 0.0) + components.get("ice_x", 0.0))
    ammonia = components.get("nh3_ice", 0.0)
    if water <= 0.0 and ammonia <= 0.0:
        return None
    return state("water" if water >= ammonia else "ammonia", pressure_gpa, temperature_k)
