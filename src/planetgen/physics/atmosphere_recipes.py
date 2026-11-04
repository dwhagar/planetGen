# planetgen/physics/atmosphere_recipes.py

"""
The atmosphere possibility tables (GEN.478, docs/design/composition/atmosphere-composition-parts/09-weighted-composition-tables.md)
========================================================================================================================

The recipes of every atmosphere, cloud, trace gas, ground radioisotope and host-star UV regime as gated rows in the format of
the possibility table (`generation/possibility_table.py`, GEN.472): 16 tables a to p and 104 rows, each with the kinds of body
and conditions of its gate, a weight that may depend on the mass, temperature, pressure, redox, insolation, the host star and the
formation distance over the snow line, the source of the weight with its verification label, and a source state. The words
are part 09 of the research; the machine columns are `scripts/atmosphere_recipes_data.py`; the file is
`physics/data/atmosphere_recipes.json`, built by `scripts/build_atmosphere_recipes.py` and checked in.

This item holds the tables and the draw of their rows and changes nothing the generator stores: no generator draws from them
yet (GEN.411, GEN.413, GEN.414, GEN.451, GEN.454 to GEN.456 do), so it takes no generator epoch of its own and rides the one
bump of the 8.4 stored-value chain (Boss, 6:08 PM EDT 2026-10-18, fork 1).

Kinds of body (GEN.840)
-----------------------

A row's gate names the kinds of body it supports (`BODY_KINDS`: terrestrial, giant, brown-dwarf, asteroid, comet) and the stored
conditions that must also hold; it never names a class letter. The class is derived from the finished body (`planet_classes.classify`),
so a rule such as "a world with an air" is the condition `pressure_bar >= 1e-4` and "a hot rocky world" is its zone and its dayside
temperature. A planet or moon is a `terrestrial` or a `giant` by its stored body type.

Rows and kinds
--------------

- `option`: one of the rows of a *group* that compete for 100 percent (a5 to a7, b1 to b4, ...). The rows that pass the gate
  are weighted, normalised and one is drawn; if every weight is zero the group's stored `fallback` row wins (the airless,
  bare or absent row; `"none"` is the explicit outcome that no row applies). In a group whose fallback is `"none"` (j, k) the
  weights are probabilities of occurrence, not shares: they are divided by their sum only when it exceeds 1, and what is left of
  1 is the explicit none (a lone row of weight 0.1 is drawn one time in ten). A row drawn under a parent row is gated on
  `drawn.<group>`: the draw of a group stores `drawn.<group>` in the subject, so the level-2 rows (c, e5, e6, h4, m5 ...) pass
  only after their parent won and the joint probability is parent times conditional.
- `binary`: present with probability `weight`, absent with the rest (the trace gases, the cloud classes by temperature).
- `param`: a continuous quantity as a parametric draw (`dist`: normal in the value or in its log10, uniform, log-uniform or a
  mixture) and as a binned table (`bins`): `draw_value`, `bins`. A row may set a stored quantity (`sets`) for the tables after it.
- `bins`, `derived`, `flag`, `windows`, `choice`: a binned view of another row, quantities computed from a formula, a
  0 or 1 gate, the cloud windows a temperature profile crosses, and a weighted pick among named alternatives.

A weight that is a DESIGN bracket draws a seeded constant inside it, once per table draw and in the file's order (log-uniform,
or uniform where the bracket starts at 0), and `p()` and the quick views use the geometric midpoint when the subject carries no
drawn constant. Nothing is fixed: a constant is the seeded draw or the subject's own `const:<scope>:<name>`.

Source states (Boss's rule of 3:02 AM EDT, composition-model.md section 4a item 10): a row whose weight source is read and
computed is `verified`; one that also leans on a second-hand, snippet, recalled, not-reached or design part is `partly
verified` with that as its limit; a row with no sourced weight is `open` and names its research item (R1 to R15 = GEN.480
to GEN.494); a row built from a measured Solar System body is a `stand-in` (a2, k7). The tables register in
`possibility_table.REGISTRY` as `atmosphere-a` to `atmosphere-p` so the report lists every open row with its item. Their
weights for the sourced rows are the same expressions; the open rows carry `open` there, so the draw of this module
(`Recipes.draw_group`, which gives an open row its seeded constant) is the one to use for them.

The corrections of the later research are applied (the part's section 5 and tunable-gaps.md): the C/O recipe of part 05 that
put 13.9 percent of giants above 0.8 is the explicit core and tail of f7, the lava-world boundaries are the text edges of
[R3059] (the difference to the caption list is 63 K and not 130 K), the shoreline of d is the Berta-Thompson logistic and
not the K window fitted to 13 Solar System bodies, a3 and a4 follow the 22 polluted white dwarfs, f7 has the planet tail and
the core centre of the data, and p4 depends on the stellar mass.

Use:

    RECIPES.draw_group("d", Subject("terrestrial", {...}))            # a GroupDraw: the row, its probability, all probabilities
    RECIPES.draw_value("a1", subject)                       # a continuous quantity
    RECIPES.bins("f7", subject)                             # [(low, high, probability)]
    RECIPES.draw_all(subject)                               # every table in the draw order
    RECIPES.entry("c1").state, .items, .source              # the row
"""

import ast
import json
import math
import os
import statistics
from dataclasses import dataclass, field

from planetgen import tuning
from planetgen.generation import possibility_table as pt
from planetgen.generation import tuning_harness
from planetgen.physics import comet_species, elements, envelope_escape, magma_ocean_redox, saturation_time
from planetgen.physics.elements import DataError
from planetgen.util import draw

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "atmosphere_recipes.json")
"""str: The data file."""

DRAW_ORDER = tuple("abdcefghijklmnop")
"""tuple: The order the tables are drawn in (part 09 section 1): the host star first (not a table), then a, b, d, c, e to l,
then m, n, o and p, which read earlier draws."""

BODY_KINDS = ("terrestrial", "giant", "brown-dwarf", "asteroid", "comet")
"""tuple: The kinds of body a row's gate may name (GEN.840). A planet or moon is `terrestrial` or `giant` by its stored body type
(`body_type` "t" or "g"), the class letter being derived from the finished body and so never a gate; the rest are the object kinds."""

KINDS = ("option", "binary", "param", "bins", "derived", "flag", "windows", "choice")

INPUTS = {
    "mass_earth": "mass of the body (Earth masses)",
    "mass_jupiter": "mass of the body (Jupiter masses)",
    "core_radius_earth": "radius of the rocky core (Earth radii)",
    "radius_km": "radius of the body (km)",
    "v_esc_kms": "escape speed (km/s)",
    "rotation_percentile": "initial rotation percentile of the host star (0.02 to 0.98), drawn by row p6",
    "teq_k": "equilibrium temperature (K)",
    "ts_k": "surface temperature (K)",
    "photosphere_bar": "pressure of the layer a giant planet's regime is read at (bar); about 0.1 for the dayside emission",
    "tday_k": "dayside temperature (K)",
    "tp_k": "dayside surface temperature of a lava world (K)",
    "t_layer_k": "temperature in the 0.1 to 1 bar layer (K)",
    "teff_k": "effective temperature of a brown dwarf (K)",
    "pressure_bar": "surface pressure (bar)",
    "redox_diw": "log10 fO2 above the iron-wuestite buffer of the outgassing stage (set by a1, a3 or a4)",
    "redox_diw_later": "the same for the later upper mantle (set by a2)",
    "insolation": "bolometric flux relative to the Earth's",
    "luminosity_solar": "luminosity of the host (solar)",
    "xuv_index": "insolation times the star's cumulative XUV ratio to the Sun's",
    "snow_ratio": "formation distance over the snow line, x = a_form / a_snow",
    "zone": "orbit zone: h hot, e ecosphere, c cold, r rogue",
    "host_type": "F, G, K, M, WD or PSR",
    "host_mass_solar": "mass of the host (solar)",
    "age_gy": "age (Gyr)",
    "is_moon": "True for a moon",
    "long_lived_geology": "True for a tectonic or long-lived body",
    "ocean_present": "True when an ocean can condense",
    "water_lost": "True when the water was lost",
    "runaway_ratio": "insolation over the runaway-greenhouse limit of the host",
    "perihelion_au": "perihelion distance of a comet (AU)",
    "family": "asteroid-field family (carbonaceous, stony, metallic, icy, basaltic, mixed, dust, collisional)",
    "o2_pal": "oxygen partial pressure as a fraction of the present atmospheric level",
    "co2_rich": "True for a CO2-rich atmosphere",
    "h2_poor": "True for an H2-poor atmosphere",
    "life_chemical": "True when the life stage is chemical or later",
    "sep_reducing": "True for weakly reducing air driven by stellar energetic particles",
    "lightning": "True when lightning acts",
    "gas_has_co_h2s": "True when the gas holds CO and H2S",
    "has_solid_ground": "True for a solid ground surface",
    "h2o_present": "True when the body has water",
    "carbon_bearing": "True for carbon-bearing air",
    "u_star": "friction velocity of the wind (m/s)",
    "u_star_threshold": "threshold friction velocity for dust (m/s)",
    "sulfur_flux_earth": "outgassed sulfur flux over the Earth's volcanic flux",
    "h2_air": "True for an H2-dominated atmosphere",
    "thin_disc_fgk": "True for a thin-disc F, G or K host",
    "th_u_f": "thorium-to-uranium enrichment factor of the host's birth gas (physics.rprocess_enrichment, at least 1)",
    "planetesimals_early": "True when the planetesimals formed within a few Myr",
    "grain_um": "cloud grain radius (micrometres)",
    "profile_t_min_k": "coldest temperature of the atmospheric profile (K)",
    "profile_t_max_k": "hottest temperature of the atmospheric profile (K)",
    "metallicity": "metallicity relative to solar (set by f5)",
    "c_to_o": "carbon-to-oxygen ratio (set by f7)",
    "n_over_c": "the host's nitrogen-to-carbon number ratio over the proto-solar one (1 for a solar composition; stellar_abundances.number_ratio)",
    "ar_over_c": "the host's argon-to-carbon number ratio over the proto-solar one",
    "s_over_c": "the host's sulfur-to-carbon number ratio over the proto-solar one",
}
"""dict: Every stored quantity a gate or a weight may read, with its meaning and unit. `body_class` is the class of the subject."""

FUNCTIONS = ("phi", "window", "exp", "log", "log10", "sqrt", "min", "max", "abs", "clamp", "loglog", "envelope_lost_fraction", "saturation_time_myr", "magma_ocean_diw", "magma_ocean_diw_sigma", "later_mantle_diw", "comet_species_median", "comet_species_sigma", "w", "p", "pdist")
SHAPES = ("log", "uniform")
DIST_KINDS = ("normal", "uniform", "loguniform", "mixture")

_NODES = (ast.Expression, ast.BinOp, ast.UnaryOp, ast.BoolOp, ast.Compare, ast.IfExp, ast.Call, ast.Name, ast.Constant,
          ast.Tuple, ast.Load, ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow, ast.USub, ast.UAdd, ast.Not, ast.And, ast.Or,
          ast.Lt, ast.LtE, ast.Gt, ast.GtE, ast.Eq, ast.NotEq, ast.In, ast.NotIn)


# ---------------------------------------------------------------------------------------------------------------------
# Expressions
# ---------------------------------------------------------------------------------------------------------------------

def phi(x):
    """The standard normal CDF."""
    return 0.5 * math.erfc(-x / math.sqrt(2.0))


def window(x, lo, hi, s):
    """A soft window: Phi((hi - x)/s) - Phi((lo - x)/s); a missing high end counts as 1 and a missing low end as 0."""
    upper = 1.0 if hi is None else phi((hi - x) / s)
    lower = 0.0 if lo is None else phi((lo - x) / s)
    return upper - lower


def clamp(x, lo, hi):
    return max(lo, min(hi, x))


def loglog(x, *anchors):
    """log10 y linear in log10 x through the anchor points (x1, y1, x2, y2, ...), held constant beyond the ends."""
    points = sorted(zip(anchors[0::2], anchors[1::2]))
    if x <= points[0][0]:
        return points[0][1]
    if x >= points[-1][0]:
        return points[-1][1]
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if x0 <= x <= x1:
            t = (math.log10(x) - math.log10(x0)) / (math.log10(x1) - math.log10(x0))
            return 10.0 ** (math.log10(y0) + t * (math.log10(y1) - math.log10(y0)))
    raise ValueError("no anchor interval")


_PLAIN = {"phi": phi, "window": window, "exp": math.exp, "log": math.log, "log10": math.log10, "sqrt": math.sqrt,
          "min": min, "max": max, "abs": abs, "clamp": clamp, "loglog": loglog, "envelope_lost_fraction": envelope_escape.lost_fraction,
          "saturation_time_myr": saturation_time.saturation_time_myr,
          "magma_ocean_diw": magma_ocean_redox.surface_diw, "magma_ocean_diw_sigma": magma_ocean_redox.diw_sigma,
          "later_mantle_diw": magma_ocean_redox.later_mantle_diw,
          "comet_species_median": comet_species.median_over_h2o, "comet_species_sigma": comet_species.sigma_dex}


class Expression:
    """
    A parsed weight or value expression.

    Attributes:
        text (str): The source.
        names (frozenset): The names it reads.
        references (frozenset): The row ids named in `w`, `p` and `pdist`.
    """

    def __init__(self, text, owner):
        self.text = text
        self.owner = owner
        try:
            self.tree = ast.parse(str(text), mode="eval")
        except SyntaxError as error:
            raise DataError(f"{owner}: the expression {text!r} cannot be read ({error.msg})") from None
        names, references = set(), set()
        for node in ast.walk(self.tree):
            if not isinstance(node, _NODES):
                raise DataError(f"{owner}: the expression {text!r} uses {type(node).__name__}, which is not allowed")
            if isinstance(node, ast.Call):
                if not isinstance(node.func, ast.Name) or node.func.id not in FUNCTIONS or node.keywords:
                    raise DataError(f"{owner}: the expression {text!r} calls something other than {', '.join(FUNCTIONS)}")
                if node.func.id in ("w", "p", "pdist"):
                    first = node.args[0] if node.args else None
                    if not (isinstance(first, ast.Constant) and isinstance(first.value, str)):
                        raise DataError(f"{owner}: {node.func.id}() in {text!r} takes a row id as its first argument")
                    references.add(first.value)
            elif isinstance(node, ast.Name) and node.id not in FUNCTIONS:
                names.add(node.id)
            elif isinstance(node, ast.Constant) and not isinstance(node.value, (int, float, str, bool, type(None))):
                raise DataError(f"{owner}: the expression {text!r} holds a constant of an unsupported type")
        self.names = frozenset(names)
        self.references = frozenset(references)

    def __call__(self, env):
        """Evaluates it; `env` has `get(name)` (lazy) and `call(name, args)` for w, p and pdist."""
        return self._eval(self.tree.body, env)

    def _eval(self, node, env):
        if isinstance(node, ast.Constant):
            return node.value
        if isinstance(node, ast.Name):
            return env.get(node.id)
        if isinstance(node, ast.Tuple):
            return tuple(self._eval(e, env) for e in node.elts)
        if isinstance(node, ast.UnaryOp):
            value = self._eval(node.operand, env)
            return {ast.USub: lambda v: -v, ast.UAdd: lambda v: +v, ast.Not: lambda v: not v}[type(node.op)](value)
        if isinstance(node, ast.BinOp):
            left, right = self._eval(node.left, env), self._eval(node.right, env)
            kind = type(node.op)
            if kind is ast.Add:
                return left + right
            if kind is ast.Sub:
                return left - right
            if kind is ast.Mult:
                return left * right
            if kind is ast.Div:
                return left / right
            try:
                return left ** right
            except OverflowError:
                return math.inf
        if isinstance(node, ast.BoolOp):
            result = None
            for operand in node.values:
                result = self._eval(operand, env)
                if isinstance(node.op, ast.And) and not result:
                    return result
                if isinstance(node.op, ast.Or) and result:
                    return result
            return result
        if isinstance(node, ast.Compare):
            left = self._eval(node.left, env)
            for operator, comparator in zip(node.ops, node.comparators):
                right = self._eval(comparator, env)
                ok = {ast.Lt: lambda a, b: a < b, ast.LtE: lambda a, b: a <= b, ast.Gt: lambda a, b: a > b,
                      ast.GtE: lambda a, b: a >= b, ast.Eq: lambda a, b: a == b, ast.NotEq: lambda a, b: a != b,
                      ast.In: lambda a, b: a in b, ast.NotIn: lambda a, b: a not in b}[type(operator)](left, right)
                if not ok:
                    return False
                left = right
            return True
        if isinstance(node, ast.IfExp):
            return self._eval(node.body if self._eval(node.test, env) else node.orelse, env)
        if isinstance(node, ast.Call):
            args = [self._eval(a, env) for a in node.args]
            name = node.func.id
            if name in _PLAIN:
                return _PLAIN[name](*args)
            return env.call(name, args)
        raise DataError(f"{self.owner}: cannot evaluate {type(node).__name__}")


class _Env:
    """The names an expression reads at one subject."""

    def __init__(self, recipes, entry, subject, lets=None, value=None, stack=()):
        self.recipes = recipes
        self.entry = entry
        self.subject = subject
        self.lets = lets if lets is not None else {}
        self.value = value
        self.stack = stack

    def get(self, name):
        if name in self.lets:
            return self.lets[name]
        if name == "value" and self.value is not None:
            return self.value
        if name == "body_class":
            return self.subject.body_class
        scope = self.entry.constants
        if name in scope:
            return self.recipes.constant(self.entry.id, name, scope[name], self.subject)
        group = self.recipes.groups.get(self.entry.group) if self.entry.group else None
        if group is not None and name in group.constants:
            return self.recipes.constant(group.name, name, group.constants[name], self.subject)
        if name in self.recipes.physical:
            return self.recipes.physical[name]
        return self.subject[name]

    def call(self, name, args):
        if name == "w":
            return self.recipes.raw_weight(self.recipes.entry(args[0]), self.subject, stack=self.stack)
        if name == "p":
            return self.recipes.probability_of(args[0], self.subject, stack=self.stack)
        if name == "pdist":
            return self.recipes.probability_in(args[0], self.subject, args[1], args[2], stack=self.stack)
        raise DataError(f"unknown function {name}")


# ---------------------------------------------------------------------------------------------------------------------
# The rows
# ---------------------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class Entry:
    """
    One row of part 09.

    Attributes:
        id (str): `a1` to `p6`.
        category (str): The table letter.
        kind (str): One of `KINDS`.
        group (str or None): The group an `option` competes in.
        option (str): What the row is, in the part's words.
        gate_text, weight_text (str): The part's words for the gate and the weight.
        source (str): Where the weight comes from, with its verification labels.
        labels (tuple): The labels in the source.
        status (str): SOURCED, DERIVED or DESIGN as the part has it.
        state (str): The source state in the loader's words (`possibility_table.STATES`).
        limit (str): For a partly verified row, its limit.
        replacement (str): For a stand-in, the builder item (GEN.509 to GEN.525, GEN.595) that replaces it.
        research (tuple): For an open row, its research ids (R1 to R15).
        items (tuple): The TODO items that replace it (GEN.480 to GEN.494 for an open row).
        classes (tuple): The kinds of body the row supports (`BODY_KINDS`), or `("*",)`.
        conditions (tuple): `Condition`s.
        weight (Expression or None): The weight (option, binary, flag).
        lets (dict): Named intermediate expressions, evaluated in order before the weight.
        constants (dict): The design constants of this row: name -> (low, high, shape).
        notes (str): The part's notes with the corrections.
        tunable (bool): The row's knobs are tuned to a published statistic: `physics/data/tuning_targets.json` must hold a
            target for it (GEN.503).
        spec (dict): The rest of the row as the file has it (`dist`, `dists`, `bins`, `values`, `windows`, `options`, `table`, `of`, `sets`).
    """

    id: str
    category: str
    kind: str
    group: str | None
    option: str
    gate_text: str
    weight_text: str
    source: str
    labels: tuple
    status: str
    state: str
    limit: str
    replacement: str
    research: tuple
    items: tuple
    classes: tuple
    conditions: tuple
    weight: object
    lets: dict
    constants: dict
    notes: str
    tunable: bool = False
    spec: dict = field(default_factory=dict)

    def supports(self, subject):
        if pt.ANY_CLASS not in self.classes and subject.body_class not in self.classes:
            return False
        for condition in self.conditions:
            if condition.quantity.startswith("drawn.") and subject.get(condition.quantity) is None:
                return False
            if not condition.holds(subject):
                return False
        return True


@dataclass(frozen=True)
class Group:
    name: str
    rows: tuple
    fallback: str | None
    constants: dict


@dataclass(frozen=True)
class GroupDraw:
    """
    A drawn group.

    Attributes:
        group (str): The group.
        row (Entry or None): The row drawn; `None` for the explicit "none applies" fallback.
        probability (float): Its probability among the rows that passed (1 for a fallback).
        probabilities (tuple): `((id, probability), ...)` of every row that passed.
        constants (dict): The design constants drawn for this draw, `{(scope, name): value}`.
        fallback (bool): True when every weight was zero and the group's fallback won.
    """

    group: str
    row: object
    probability: float
    probabilities: tuple
    constants: dict
    fallback: bool = False

    @property
    def name(self):
        return self.row.id if self.row is not None else "none"


@dataclass
class RecipeDraw:
    """What `draw_all` drew: groups, values, binaries, windows and choices, and the joint probability of every group."""

    groups: dict = field(default_factory=dict)
    values: dict = field(default_factory=dict)
    binaries: dict = field(default_factory=dict)
    windows: dict = field(default_factory=dict)
    choices: dict = field(default_factory=dict)
    flags: dict = field(default_factory=dict)
    joint: dict = field(default_factory=dict)
    subject: object = None


class _PtWeight:
    """A weight function for the registered possibility table; its qualified name carries the expression so the digest sees it."""

    def __init__(self, recipes, entry):
        self.recipes = recipes
        self.entry = entry
        self.__module__ = __name__
        self.__qualname__ = f"{entry.id}={entry.weight.text if entry.weight is not None else ''}"

    def __call__(self, subject):
        return self.recipes.raw_weight(self.entry, subject, gated=False)


class Recipes:
    """The 104 rows, their groups and the draws."""

    def __init__(self, document, register=True):
        self.document = document
        self.title = document["title"]
        self.physical = dict(document["constants"])
        self.research = dict(document["research"])
        self.categories = dict(document["categories"])
        self.groups = {name: Group(name, tuple(spec["rows"]), spec["fallback"],
                                   {k: self._constant_spec(f"group {name}", k, v)
                                    for k, v in spec.get("constants", {}).items()})
                       for name, spec in document["groups"].items()}
        self.entries = []
        self._by_id = {}
        for raw in document["entries"]:
            entry = self._entry(raw)
            if entry.id in self._by_id:
                raise DataError(f"row {entry.id} appears twice")
            self.entries.append(entry)
            self._by_id[entry.id] = entry
        self._check()
        self.tables = self._register() if register else self._build_tables()

    # -- loading ---------------------------------------------------------------------------------------------------

    def _entry(self, raw):
        row_id = raw.get("id")
        if not row_id or not isinstance(row_id, str):
            raise DataError("a row has no id")
        for key in ("category", "kind", "option", "source", "status", "state", "classes"):
            if not raw.get(key):
                raise DataError(f"row {row_id}: no {key}")
        if raw["kind"] not in KINDS:
            raise DataError(f"row {row_id}: kind {raw['kind']!r} is none of {', '.join(KINDS)}")
        if raw["category"] not in self.categories:
            raise DataError(f"row {row_id}: category {raw['category']!r} is unknown")
        if not str(raw["source"]).strip():
            raise DataError(f"row {row_id} has no source")
        if raw["state"] not in pt.STATES:
            raise DataError(f"row {row_id}: state {raw['state']!r} is none of {', '.join(pt.STATES)}")
        if raw["state"] == pt.PARTLY_VERIFIED and not raw.get("limit"):
            raise DataError(f"row {row_id} is partly verified and names no limit")
        if raw["state"] == pt.STAND_IN and not raw.get("replacement"):
            raise DataError(f"row {row_id} is a stand-in and names no replacement item")
        if raw["state"] == pt.OPEN_STATE:
            if not raw.get("items") or any(item not in self.research.values() for item in raw["items"]):
                raise DataError(f"row {row_id} is open and names no research item of the part (R1 to R15)")
        known = set(BODY_KINDS) | {pt.ANY_CLASS}
        unknown = sorted(set(raw["classes"]) - known)
        if unknown:
            raise DataError(f"row {row_id}: unknown kinds {', '.join(unknown)}")
        conditions = []
        for item in raw.get("conditions", []):
            extra = set(item) - {"quantity", "minimum", "maximum", "equals", "one_of"}
            if extra or "quantity" not in item:
                raise DataError(f"row {row_id}: a condition with {sorted(extra) or 'no quantity'}")
            quantity = item["quantity"]
            if quantity not in INPUTS and not quantity.startswith("drawn."):
                raise DataError(f"row {row_id}: the condition reads {quantity!r}, which is no stored quantity")
            kwargs = dict(item)
            if "one_of" in kwargs:
                kwargs["one_of"] = tuple(kwargs["one_of"])
            conditions.append(pt.Condition(**kwargs))
        constants = {}
        for name, spec in raw.get("constants", {}).items():
            constants[name] = self._constant_spec(row_id, name, spec)
        lets = {name: Expression(text, f"row {row_id} let {name}") for name, text in raw.get("lets", {}).items()}
        weight = raw.get("weight")
        expression = None
        if weight is not None:
            expression = Expression(str(weight), f"row {row_id}")
        if raw["kind"] in ("option", "binary", "flag") and expression is None:
            raise DataError(f"row {row_id} is an {raw['kind']} and has no weight")
        spec = {k: raw[k] for k in ("dist", "dists", "bins", "values", "windows", "options", "table", "table_of", "of",
                                    "sets") if k in raw}
        return Entry(
            id=row_id, category=raw["category"], kind=raw["kind"], group=raw.get("group"), option=raw["option"],
            gate_text=raw.get("gate_text", ""), weight_text=raw.get("weight_text", ""), source=raw["source"],
            labels=tuple(raw.get("labels", [])), status=raw["status"], state=raw["state"], limit=raw.get("limit", ""),
            replacement=raw.get("replacement", ""), research=tuple(raw.get("research", [])),
            items=tuple(raw.get("items", [])), classes=tuple(raw["classes"]), conditions=tuple(conditions),
            weight=expression, lets=lets, constants=constants, notes=raw.get("notes", ""),
            tunable=bool(raw.get("tunable", False)), spec=spec)

    @staticmethod
    def _constant_spec(owner, name, spec):
        if not (isinstance(spec, (list, tuple)) and len(spec) == 3):
            raise DataError(f"{owner}: the constant {name} is not [low, high, shape]")
        low, high, shape = spec
        if shape not in SHAPES or not (isinstance(low, (int, float)) and isinstance(high, (int, float)) and low <= high):
            raise DataError(f"{owner}: the constant {name} is not a bracket low <= high with a shape of {', '.join(SHAPES)}")
        if shape == "log" and low <= 0:
            raise DataError(f"{owner}: the log-uniform constant {name} starts at {low}")
        return (float(low), float(high), shape)

    def _check(self):
        for name, group in self.groups.items():
            for row_id in group.rows:
                if row_id not in self._by_id or self._by_id[row_id].group != name:
                    raise DataError(f"group {name}: the row {row_id} is missing or belongs elsewhere")
            if group.fallback not in (None, "none") and group.fallback not in group.rows:
                raise DataError(f"group {name}: the fallback {group.fallback} is not one of its rows")
        for entry in self.entries:
            if entry.kind == "option" and (entry.group not in self.groups or entry.id not in self.groups[entry.group].rows):
                raise DataError(f"row {entry.id} is an option of no group")
            self._check_names(entry)
            self._check_spec(entry)

    def _check_names(self, entry):
        allowed = set(INPUTS) | set(self.physical) | {"body_class", "value"} | set(entry.constants) | set(entry.lets)
        if entry.group:
            allowed |= set(self.groups[entry.group].constants)
        expressions = list(entry.lets.values())
        if entry.weight is not None:
            expressions.append(entry.weight)
        spec = entry.spec
        for key in ("values",):
            for name, text in spec.get(key, {}).items():
                expressions.append(Expression(text, f"row {entry.id} value {name}"))
        for expression in self._dist_expressions(entry):
            expressions.append(expression)
        for expression in expressions:
            unknown = sorted(expression.names - allowed)
            if unknown:
                raise DataError(f"row {entry.id}: {expression.text!r} reads {', '.join(unknown)}, which is no stored "
                                f"quantity, constant or intermediate")
            for reference in expression.references:
                if reference not in self._by_id:
                    raise DataError(f"row {entry.id}: {expression.text!r} names the row {reference}, which does not exist")

    def _dist_expressions(self, entry):
        out = []

        def walk(dist, label):
            for key in ("mu", "sigma", "lo", "hi", "weight"):
                if key in dist and isinstance(dist[key], str):
                    out.append(Expression(dist[key], f"row {entry.id} {label} {key}"))
            for bound in dist.get("range", []) or []:
                if isinstance(bound, str):
                    out.append(Expression(bound, f"row {entry.id} {label} range"))
            for sub in ("core", "tail"):
                if sub in dist:
                    walk(dist[sub], f"{label} {sub}")

        spec = entry.spec
        if "dist" in spec:
            walk(spec["dist"], "dist")
        for name, dist in spec.get("dists", {}).items():
            walk(dist, name)
        for edge in spec.get("bins", []):
            if isinstance(edge, str):
                out.append(Expression(edge, f"row {entry.id} bins"))
        return out

    def _check_spec(self, entry):
        spec, kind = entry.spec, entry.kind
        if kind == "param" and not (("dist" in spec) ^ ("dists" in spec)):
            raise DataError(f"row {entry.id} is a param and needs exactly one of dist and dists")
        for dist in ([spec["dist"]] if "dist" in spec else []) + list(spec.get("dists", {}).values()):
            self._check_dist(entry.id, dist)
        if kind == "binary" and "dist" in spec:
            self._check_dist(entry.id, spec["dist"])
        if kind == "windows":
            if not spec.get("windows"):
                raise DataError(f"row {entry.id} lists no windows")
            for name, low, high in spec["windows"]:
                if not (isinstance(low, (int, float)) and isinstance(high, (int, float)) and low < high):
                    raise DataError(f"row {entry.id}: the window {name} is not low < high")
        if kind == "choice":
            options = spec.get("options") or {}
            if not options or not all(o.get("weight", 0) > 0 for o in options.values()):
                raise DataError(f"row {entry.id}: a choice needs options with a positive weight")
        if kind == "bins":
            if spec.get("of") not in self._by_id or self._by_id[spec["of"]].kind != "param":
                raise DataError(f"row {entry.id}: a bins row names a param row it is the bins of")
        if kind == "derived" and not (spec.get("values") or spec.get("table")):
            raise DataError(f"row {entry.id} is derived and has no values")
        if "table" in spec:
            group = self.groups.get(spec.get("table_of"))
            if group is None or set(spec["table"]) != set(group.rows):
                raise DataError(f"row {entry.id}: the table must have a value for every row of group {spec.get('table_of')}")
        if entry.spec.get("sets") and entry.spec["sets"] not in INPUTS:
            raise DataError(f"row {entry.id} sets {entry.spec['sets']!r}, which is no stored quantity")
        bins = spec.get("bins")
        if bins is not None and kind == "param":
            numeric = [b for b in bins if isinstance(b, (int, float))]
            if len(numeric) == len(bins) and sorted(numeric) != list(numeric):
                raise DataError(f"row {entry.id}: the bin edges are not ascending")

    @staticmethod
    def _check_dist(owner, dist):
        kind = dist.get("kind")
        if kind not in DIST_KINDS:
            raise DataError(f"row {owner}: distribution kind {kind!r} is none of {', '.join(DIST_KINDS)}")
        need = {"normal": ("mu", "sigma"), "uniform": ("lo", "hi"), "loguniform": ("lo", "hi"),
                "mixture": ("weight", "core", "tail")}[kind]
        missing = [k for k in need if k not in dist]
        if missing:
            raise DataError(f"row {owner}: the {kind} distribution lacks {', '.join(missing)}")
        if kind == "mixture":
            for sub in ("core", "tail"):
                Recipes._check_dist(owner, dist[sub])

    def _build_tables(self):
        tables = []
        for letter in sorted(self.categories):
            rows = [self._pt_row(e) for e in self.entries if e.category == letter]
            tables.append(pt.Table(f"atmosphere-{letter}", self.document["version"], rows,
                                   description=f"{self.categories[letter]} (part 09 category {letter})"))
        return tuple(tables)

    def _register(self):
        tables = self._build_tables()
        tuning_harness.TARGETS.require(tables)
        return tuple(pt.REGISTRY.register(table) for table in tables)

    def _pt_row(self, entry):
        if entry.state == pt.OPEN_STATE:
            weight, source, notes = pt.OPEN, pt.OPEN, f"{entry.option}. Source: {entry.source}. {entry.notes}".strip()
        else:
            weight = _PtWeight(self, entry) if entry.kind in ("option", "binary", "flag") else 1.0
            source, notes = entry.source, f"{entry.option}. {entry.notes}".strip()
        return pt.Row(name=entry.id, classes=entry.classes, weight=weight, source=source, state=entry.state,
                      conditions=entry.conditions, limit=entry.limit, replacement=entry.replacement, notes=notes,
                      tunable=entry.tunable)

    # -- access ----------------------------------------------------------------------------------------------------

    def entry(self, row_id):
        try:
            return self._by_id[row_id]
        except KeyError:
            raise KeyError(f"no atmosphere row {row_id!r}") from None

    def __len__(self):
        return len(self.entries)

    def in_category(self, letter):
        return [e for e in self.entries if e.category == letter]

    def state_counts(self):
        counts = {state: 0 for state in pt.STATES}
        for entry in self.entries:
            counts[entry.state] += 1
        return counts

    def open_rows(self):
        """`[(row, research ids, items)]` of the open rows: what each waits for."""
        return [(e, e.research, e.items) for e in self.entries if e.state == pt.OPEN_STATE]

    # -- constants -------------------------------------------------------------------------------------------------

    @staticmethod
    def _midpoint(bracket):
        low, high, shape = bracket
        if shape == "log" and low > 0:
            return math.sqrt(low * high)
        return 0.5 * (low + high)

    def constant(self, scope, name, bracket, subject):
        """The design constant: the subject's own `const:<scope>:<name>`, else the bracket's midpoint (geometric for log)."""
        value = subject.get(f"const:{scope}:{name}")
        return self._midpoint(bracket) if value is None else value

    def _draw_constants(self, scope, constants):
        drawn = {}
        for name, (low, high, shape) in constants.items():
            if low == high:
                value = low
            elif shape == "log":
                value = math.exp(draw.uniform(math.log(low), math.log(high)))
            else:
                value = draw.uniform(low, high)
            drawn[(scope, name)] = value
        return drawn

    def _with_constants(self, subject, drawn):
        quantities = dict(subject.quantities)
        for (scope, name), value in drawn.items():
            quantities[f"const:{scope}:{name}"] = value
        return pt.Subject(subject.body_class, quantities)

    # -- weights ---------------------------------------------------------------------------------------------------

    def _env(self, entry, subject, stack, value=None):
        env = _Env(self, entry, subject, stack=stack, value=value)
        for name, expression in entry.lets.items():
            env.lets[name] = self._safe(entry, expression, env)
        return env

    def _safe(self, entry, expression, env):
        try:
            return expression(env)
        except (ZeroDivisionError, ValueError, OverflowError, TypeError) as error:
            raise pt.TableError(f"row {entry.id}: {expression.text!r} cannot be evaluated here ({error})") from None

    def raw_weight(self, entry, subject, gated=True, stack=()):
        """The weight of `entry` at `subject`: 0 where the gate fails (unless `gated` is False), else the evaluated expression."""
        if entry.id in stack:
            raise pt.TableError(f"row {entry.id}: weights refer to each other in a circle ({' -> '.join(stack)})")
        if gated and not entry.supports(subject):
            return 0.0
        if entry.weight is None:
            return 1.0
        stack = stack + (entry.id,)
        env = self._env(entry, subject, stack)
        value = self._safe(entry, entry.weight, env)
        if isinstance(value, bool):
            value = float(value)
        if not isinstance(value, (int, float)) or not math.isfinite(value) or value < -1e-12:
            raise pt.TableError(f"row {entry.id}: the weight evaluated to {value!r}, not a finite number of 0 or more")
        value = max(0.0, float(value))
        if entry.kind in ("binary", "flag") and value > 1.0 + 1e-12:
            raise pt.TableError(f"row {entry.id}: a probability of {value} is above 1")
        return value

    def _class_in_group(self, name, subject):
        return any(pt.ANY_CLASS in e.classes or subject.body_class in e.classes
                   for e in (self.entry(r) for r in self.groups[name].rows))

    def group_probabilities(self, name, subject, stack=()):
        """
        `[(row, probability)]` of the rows of group `name` that pass `subject`, in the file's order, adding up to 1 (in a group whose
        fallback is `"none"` the weights are absolute and the shortfall to 1 is a `(None, p)` entry). Where no row
        has a weight above zero (or the conditions of every row fail) the group's fallback row takes the whole (`None` for the
        explicit "none applies" outcome); `[]` where the class is not one the group supports at all, or where it has no fallback
        and nothing has a weight.
        """
        group = self.groups[name]
        if not self._class_in_group(name, subject):
            return []
        rows = [self.entry(r) for r in group.rows]
        passing = [r for r in rows if r.supports(subject)]
        weights = [(r, self.raw_weight(r, subject, gated=False, stack=stack)) for r in passing]
        weights = [(r, w) for r, w in weights if w > 0.0]
        total = math.fsum(w for _r, w in weights)
        if weights and total > 0.0 and group.fallback == "none":
            # The weights of a group whose fallback is "none applies" are probabilities of occurrence: a lone row of weight 0.1 is
            # drawn one time in ten, and the rest is the explicit none (GEN.674: it used to be normalised to 1).
            scale = max(total, 1.0)
            table = [(r, w / scale) for r, w in weights]
            if total < 1.0:
                table.append((None, 1.0 - total))
            return table
        if weights and total > 0.0:
            return [(r, w / total) for r, w in weights]
        if group.fallback is None:
            return []
        return [(None if group.fallback == "none" else self.entry(group.fallback), 1.0)]

    def probability_of(self, row_id, subject, stack=()):
        """The probability row `row_id` has in its group at `subject` (0 when the group does not apply)."""
        entry = self.entry(row_id)
        if entry.kind != "option":
            raise pt.TableError(f"{row_id} is not an option of a group, so it has no probability among rows")
        for row, probability in self.group_probabilities(entry.group, subject, stack=stack):
            if row is not None and row.id == row_id:
                return probability
        return 0.0

    # -- draws -----------------------------------------------------------------------------------------------------

    def probabilities(self, group, subject):
        """The probabilities of the rows of `group` at `subject` with each design constant at the subject's value or
        the bracket's midpoint (no random draw is made): `[(id, probability)]`."""
        table = self.group_probabilities(group, subject)
        return [(row.id if row is not None else "none", p) for row, p in table]

    def draw_group(self, group, subject):
        """
        Draws one row of `group` for `subject`: the design constants of the group's rows are drawn first (seeded, in the file's
        order), the rows that pass are weighted and normalised and one is picked with the one random wrapper.

        Raises:
            pt.PossibilityError: The class is none the group's rows support, or every weight is zero and the group has no
                fallback.
        """
        spec = self.groups[group]
        drawn = self._draw_constants(group, spec.constants)
        for row_id in spec.rows:
            drawn.update(self._draw_constants(row_id, self.entry(row_id).constants))
        view = self._with_constants(subject, drawn)
        rows = [self.entry(r) for r in spec.rows]
        if not self._class_in_group(group, view):
            raise pt.PossibilityError(f"atmosphere-{rows[0].category}", view.body_class,
                                      f"no row of the group {group} supports the class")
        table = self.group_probabilities(group, view)
        if not table:
            raise pt.PossibilityError(f"atmosphere-{rows[0].category}", view.body_class,
                                      f"no row of the group {group} has a weight above zero and it has no fallback")
        fallback = len(table) == 1 and (table[0][0] is None or table[0][0].id == spec.fallback) and not any(
            r.supports(view) and self.raw_weight(r, view, gated=False) > 0.0 for r in rows)
        if fallback:
            row = table[0][0]
            return GroupDraw(group, row, 1.0, ((row.id if row else "none", 1.0),), drawn, True)
        chosen = draw.choices([r for r, _p in table], weights=[p for _r, p in table], k=1)[0]
        listed = tuple((r.id if r else "none", p) for r, p in table)
        probability = dict(listed)[chosen.id if chosen else "none"]
        return GroupDraw(group, chosen, probability, listed, drawn)

    # -- continuous rows -------------------------------------------------------------------------------------------

    def _number(self, entry, spec, env):
        if isinstance(spec, (int, float)):
            return float(spec)
        return float(self._safe(entry, Expression(spec, f"row {entry.id}"), env))

    def _dist(self, entry, dist, env):
        """The distribution with its numbers evaluated at the subject: a dict of floats."""
        kind = dist["kind"]
        out = {"kind": kind, "log10": bool(dist.get("log10", False))}
        if kind == "normal":
            out["mu"], out["sigma"] = self._number(entry, dist["mu"], env), self._number(entry, dist["sigma"], env)
            if not out["sigma"] > 0:
                raise pt.TableError(f"row {entry.id}: the standard deviation {out['sigma']} is not positive")
            rng = dist.get("range")
            out["range"] = tuple(self._number(entry, b, env) for b in rng) if rng else None
        elif kind in ("uniform", "loguniform"):
            out["lo"], out["hi"] = self._number(entry, dist["lo"], env), self._number(entry, dist["hi"], env)
            if not out["lo"] < out["hi"] or (kind == "loguniform" and out["lo"] <= 0):
                raise pt.TableError(f"row {entry.id}: the {kind} bracket [{out['lo']}, {out['hi']}] is not valid")
        else:
            out["weight"] = self._number(entry, dist["weight"], env)
            if not 0.0 <= out["weight"] <= 1.0:
                raise pt.TableError(f"row {entry.id}: the mixture weight {out['weight']} is not in [0, 1]")
            out["core"] = self._dist(entry, dist["core"], env)
            out["tail"] = self._dist(entry, dist["tail"], env)
        return out

    @staticmethod
    def _cdf(d, v):
        kind = d["kind"]
        if kind == "mixture":
            w = d["weight"]
            return (1.0 - w) * Recipes._cdf(d["core"], v) + w * Recipes._cdf(d["tail"], v)
        if kind == "uniform":
            return min(1.0, max(0.0, (v - d["lo"]) / (d["hi"] - d["lo"])))
        if kind == "loguniform":
            if v <= d["lo"]:
                return 0.0
            return min(1.0, math.log(v / d["lo"]) / math.log(d["hi"] / d["lo"]))
        if v <= 0 and d["log10"]:
            return 0.0
        x = math.log10(v) if d["log10"] else v
        if d["range"]:
            lo, hi = d["range"]
            f = lambda t: phi((t - d["mu"]) / d["sigma"])  # noqa: E731
            low, high = f(lo), f(hi)
            if x <= lo:
                return 0.0
            if x >= hi:
                return 1.0
            return (f(x) - low) / (high - low)
        return phi((x - d["mu"]) / d["sigma"])

    @staticmethod
    def _quantile(d, u):
        kind = d["kind"]
        if kind == "uniform":
            return d["lo"] + u * (d["hi"] - d["lo"])
        if kind == "loguniform":
            return d["lo"] * (d["hi"] / d["lo"]) ** u
        norm = statistics.NormalDist()
        low, high = 0.0, 1.0
        if d["range"]:
            lo, hi = d["range"]
            low, high = phi((lo - d["mu"]) / d["sigma"]), phi((hi - d["mu"]) / d["sigma"])
        q = min(max(low + u * (high - low), 1e-16), 1.0 - 1e-16)
        x = d["mu"] + d["sigma"] * norm.inv_cdf(q)
        if d["range"]:
            x = min(max(x, d["range"][0]), d["range"][1])
        return 10.0 ** x if d["log10"] else x

    def _draw_dist(self, d):
        if d["kind"] == "mixture":
            pick = draw.random()
            part = d["tail"] if pick < d["weight"] else d["core"]
            return self._quantile(part, draw.random())
        return self._quantile(d, draw.random())

    def _dist_of(self, entry, subject):
        env = self._env(entry, subject, ())
        spec = entry.spec
        if "dists" in spec:
            return {name: self._dist(entry, d, env) for name, d in spec["dists"].items()}
        return self._dist(entry, spec["dist"], env)

    def draw_value(self, row_id, subject):
        """
        The continuous quantity of a `param` row (or the value of a `binary` row that carries one) drawn at `subject`:
        a float, or a dict of floats for a row with several (`k6`, `k7`). `None` when the row's gate does not pass.
        """
        entry = self.entry(row_id)
        if entry.kind not in ("param", "binary") or not ({"dist", "dists"} & set(entry.spec)):
            raise pt.TableError(f"row {row_id} has no distribution to draw")
        if not entry.supports(subject):
            return None
        drawn = self._draw_constants(entry.id, entry.constants)
        view = self._with_constants(subject, drawn)
        spec = self._dist_of(entry, view)
        if "dists" in entry.spec:
            return {name: self._draw_dist(d) for name, d in spec.items()}
        return self._draw_dist(spec)

    def bins(self, row_id, subject):
        """
        `[(low, high, probability)]`: the binned table of a continuous row (or of the row a `bins` row views) at `subject`.
        The ends that hold no mass are left out; the probabilities add up to 1.
        """
        entry = self.entry(row_id)
        if entry.kind == "bins":
            return self.bins(entry.spec["of"], subject) if "bins" not in entry.spec else self._bins(
                self.entry(entry.spec["of"]), entry.spec["bins"], subject)
        if "bins" not in entry.spec:
            raise pt.TableError(f"row {row_id} has no bins")
        return self._bins(entry, entry.spec["bins"], subject)

    def _bins(self, entry, edges, subject):
        drawn = self._with_constants(subject, {})
        env = self._env(entry, drawn, ())
        numbers = [self._number(entry, e, env) for e in edges]
        d = self._dist_of(entry, drawn)
        cdf = [self._cdf(d, v) for v in numbers]
        cut = [0.0] + cdf + [1.0]
        edges_ = [-math.inf] + numbers + [math.inf]
        out = []
        for i in range(len(cut) - 1):
            p = cut[i + 1] - cut[i]
            if p > 1e-12 or (0 < i < len(cut) - 2):
                out.append((edges_[i], edges_[i + 1], max(0.0, p)))
        return out

    def probability_in(self, row_id, subject, low, high, stack=()):
        """The probability that the continuous row's value lies in [low, high] (None is open)."""
        entry = self.entry(row_id)
        d = self._dist_of(entry, subject)
        lower = 0.0 if low is None else self._cdf(d, low)
        upper = 1.0 if high is None else self._cdf(d, high)
        return max(0.0, upper - lower)

    def draw_binary(self, row_id, subject):
        """`(present, value)` for a `binary` row: present with its probability (a seeded draw), the value drawn if it has one;
        `(False, None)` where its gate does not pass."""
        entry = self.entry(row_id)
        if entry.kind != "binary":
            raise pt.TableError(f"row {row_id} is not a binary row")
        if not entry.supports(subject):
            return False, None
        drawn = self._draw_constants(entry.id, entry.constants)
        view = self._with_constants(subject, drawn)
        probability = self.raw_weight(entry, view, gated=False)
        present = draw.random() < probability
        value = self._draw_dist(self._dist_of(entry, view)) if present and "dist" in entry.spec else None
        return present, value

    def probability(self, row_id, subject):
        """The presence probability of a `binary` or `flag` row at `subject` with the constants at the subject's value or the
        bracket midpoint (0 where the gate fails)."""
        entry = self.entry(row_id)
        return self.raw_weight(entry, subject)

    def crossed_windows(self, row_id, subject):
        """The names of the cloud windows of a `windows` row that the subject's temperature profile crosses."""
        entry = self.entry(row_id)
        if entry.kind != "windows":
            raise pt.TableError(f"row {row_id} is not a windows row")
        if not entry.supports(subject):
            return []
        lo_t, hi_t = subject["profile_t_min_k"], subject["profile_t_max_k"]
        return [name for name, low, high in entry.spec["windows"] if lo_t <= high and hi_t >= low]

    def derive(self, row_id, subject, value=None):
        """The named values of a `derived` row, or the extra values of another row (`value` is the drawn quantity)."""
        entry = self.entry(row_id)
        spec = entry.spec
        if "table" in spec:
            group_row = subject.get(f"drawn.{spec['table_of']}")
            if group_row not in spec["table"]:
                raise pt.TableError(f"row {row_id}: the group {spec['table_of']} has not been drawn")
            return {"p_max_bar": spec["table"][group_row]}
        if not entry.supports(subject):
            return {}
        env = self._env(entry, subject, (), value=value)
        return {name: float(self._safe(entry, Expression(text, f"row {entry.id} value {name}"), env))
                for name, text in spec.get("values", {}).items()}

    def flag(self, row_id, subject):
        """A `flag` row: True when its gate passes and its expression is 1."""
        entry = self.entry(row_id)
        if entry.kind != "flag":
            raise pt.TableError(f"row {row_id} is not a flag row")
        return self.raw_weight(entry, subject) >= 1.0

    def draw_choice(self, row_id, subject):
        """The name of the alternative a `choice` row picks (None where the gate fails)."""
        entry = self.entry(row_id)
        if not entry.supports(subject):
            return None
        names = list(entry.spec["options"])
        weights = [entry.spec["options"][n]["weight"] for n in names]
        return draw.choices(names, weights=weights, k=1)[0]

    # -- the whole draw --------------------------------------------------------------------------------------------

    def draw_all(self, subject, categories=DRAW_ORDER):
        """
        Draws the tables of `categories` for `subject` in that order, each row of a table in the file's order. A group whose
        rows none supports the subject is skipped; a group drawn stores `drawn.<group>`; a continuous row with `sets`
        stores its value under that quantity, so the tables after it read the earlier draws. Returns a `RecipeDraw`.
        """
        quantities = dict(subject.quantities)
        for name in self.groups:
            quantities.setdefault(f"drawn.{name}", None)
        current = pt.Subject(subject.body_class, quantities)
        result = RecipeDraw(subject=current)
        done = set()
        for letter in categories:
            for entry in self.in_category(letter):
                if entry.kind == "option":
                    if entry.group in done:
                        continue
                    done.add(entry.group)
                    if not any(self.entry(r).supports(current) for r in self.groups[entry.group].rows):
                        continue
                    drawn = self.draw_group(entry.group, current)
                    result.groups[entry.group] = drawn
                    current.quantities[f"drawn.{entry.group}"] = drawn.name
                    for (scope, name), value in drawn.constants.items():
                        current.quantities[f"const:{scope}:{name}"] = value
                    parent = self._parent_group(drawn.row) if drawn.row is not None else None
                    result.joint[entry.group] = drawn.probability * (result.joint.get(parent, 1.0) if parent else 1.0)
                elif entry.kind == "param":
                    value = self.draw_value(entry.id, current)
                    if value is not None:
                        result.values[entry.id] = value
                        if entry.spec.get("sets"):
                            current.quantities[entry.spec["sets"]] = value
                elif entry.kind == "binary":
                    present, value = self.draw_binary(entry.id, current)
                    if entry.supports(current):
                        result.binaries[entry.id] = (present, value)
                elif entry.kind == "windows":
                    if entry.supports(current):
                        result.windows[entry.id] = self.crossed_windows(entry.id, current)
                elif entry.kind == "choice":
                    pick = self.draw_choice(entry.id, current)
                    if pick is not None:
                        result.choices[entry.id] = pick
                elif entry.kind == "flag":
                    if entry.supports(current):
                        result.flags[entry.id] = self.flag(entry.id, current)
        return result

    def _parent_group(self, row):
        for condition in row.conditions:
            if condition.quantity.startswith("drawn."):
                return condition.quantity[len("drawn."):]
        return None


def _load():
    with open(DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    for key in ("title", "retrieved", "part", "version", "categories", "constants", "research", "groups", "entries"):
        if not document.get(key):
            raise DataError(f"atmosphere_recipes.json states no {key}")
    if not isinstance(document["version"], int) or isinstance(document["version"], bool) or document["version"] < 1:
        raise DataError("atmosphere_recipes.json: the version is a whole number from 1")
    return Recipes(document)


RECIPES = _load()
"""Recipes: The 104 rows; their 16 tables are registered in `possibility_table.REGISTRY` as `atmosphere-a` to `atmosphere-p`."""

TABLES = RECIPES.tables
