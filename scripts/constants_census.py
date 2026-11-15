"""TEST.715 (GEN.945 census by Lane 21, reports/gen945/constants_census.py, kept whole and made importable): every module-level numeric constant under
`src/planetgen`, in one class: registered (the stand-in registry or gap-filling plan names it), tagged (its own docstring or comment says stand-in, tier, [R....],
derived or read), physical (a constants or units module) or unregistered (none of these).

`python scripts/constants_census.py --write` rewrites the baseline of unregistered names (only ever shorter: the test refuses a longer one).
"""
import ast
import os
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
SOURCE = ROOT / "src" / "planetgen"
BASELINE = ROOT / "src" / "tests" / "baselines" / "unregistered_constants.txt"
GENERATOR_AREAS = ("physics", "galaxy", "generation", "population")
"""tuple: The areas the gate covers (TEST.715): the generator's own, not the engineering knobs of the web, database and queue modules. `tuning.py` is added by name."""
TAG = re.compile(r"stand-in|tier\s+T\d|\bT[1-9]\b|\[R\d{4}|\bL[0-4]\b|registry|GEN\.\d+|derived|read from|measured|published|defined|by definition|SI |IAU|CODATA", re.I)
PHYSICAL_MODULES = ("constants", "units", "unit_registry", "astropy_constants")


def _numeric(node):
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)) and not isinstance(node.value, bool):
        return True
    if isinstance(node, (ast.Tuple, ast.List)):
        return bool(node.elts) and all(_numeric(e) for e in node.elts)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        return _numeric(node.operand)
    if isinstance(node, ast.BinOp):
        return _numeric(node.left) and _numeric(node.right)
    return False


def registry_text():
    """The registry text a name is looked up in."""
    research = ROOT / "docs" / "design" / "research"
    return (research / "gap-filling-plan.md").read_text(encoding="utf-8") + (research / "stand-in-values.md").read_text(encoding="utf-8")


def census(registry=None):
    """`[{module, name, line, class}]` for every module-level numeric constant under `src/planetgen`."""
    registry = registry_text() if registry is None else registry
    words = set(re.findall(r"\w+", registry))
    dotted = set(re.findall(r"\w+\.\w+", registry))
    rows = []
    for base, _dirs, files in os.walk(SOURCE):
        for filename in sorted(files):
            if not filename.endswith(".py"):
                continue
            path = pathlib.Path(base) / filename
            text = path.read_text(encoding="utf-8")
            lines = text.split("\n")
            try:
                tree = ast.parse(text)
            except SyntaxError:
                continue
            module = ".".join(path.relative_to(ROOT / "src").with_suffix("").parts)
            for i, node in enumerate(tree.body):
                if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name) and _numeric(node.value):
                    name = node.targets[0].id
                elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name) and node.value is not None and _numeric(node.value):
                    name = node.target.id
                else:
                    continue
                if not re.fullmatch(r"_?[A-Z][A-Z0-9_]*", name):
                    continue
                end = getattr(node, "end_lineno", node.lineno)
                own = "\n".join(lines[node.lineno - 1:end])
                doc = ""
                following = tree.body[i + 1] if i + 1 < len(tree.body) else None
                if isinstance(following, ast.Expr) and isinstance(following.value, ast.Constant) and isinstance(following.value.value, str) and following.lineno == end + 1:
                    doc = following.value.value
                prev = lines[node.lineno - 2] if node.lineno >= 2 else ""
                around = own + "\n" + doc + "\n" + (prev if prev.strip().startswith("#") else "")
                if name in words or module.split(".")[-1] + "." + name in dotted:
                    kind = "registered"
                elif TAG.search(around):
                    kind = "tagged"
                elif module.endswith(PHYSICAL_MODULES):
                    kind = "physical"
                else:
                    kind = "unregistered"
                rows.append({"module": module, "name": name, "line": node.lineno, "class": kind})
    return rows


def covered(module):
    """Whether the gate covers a module: the generator's own areas and `tuning.py`."""
    parts = module.split(".")
    return parts[1] in GENERATOR_AREAS or module == "planetgen.tuning"


def unregistered(registry=None):
    """The sorted `module:name` of every uncovered-by-registry constant in the covered areas."""
    return sorted({f"{r['module']}:{r['name']}" for r in census(registry) if r["class"] == "unregistered" and covered(r["module"])})


def read_baseline():
    """The committed names, one per line."""
    return [line.strip() for line in BASELINE.read_text(encoding="utf-8").splitlines() if line.strip() and not line.startswith("#")]


HEADER = ("# TEST.715: constants of the generator's own areas (physics, galaxy, generation, population, tuning.py) that no registry row, tag or constants module covers.\n"
          "# This list may only shrink: register a name (docs/design/research/stand-in-values.md), tag it in its docstring, or move it to the constants module, and delete its line.\n")

if __name__ == "__main__":
    names = unregistered()
    if "--write" in sys.argv:
        BASELINE.parent.mkdir(parents=True, exist_ok=True)
        BASELINE.write_text(HEADER + "\n".join(names) + "\n", encoding="utf-8")
    print(len(names), "unregistered")
