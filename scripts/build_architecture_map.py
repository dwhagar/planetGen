#!/usr/bin/env python3
"""
Builds the architecture map, `docs/generated/architecture-map.html` (DOC.67).

Why this exists
===============
Boss (2026-10-17, 12:23 PM EDT): "I want a more detailed architecture HTML document that I can drill down and see
the flow through modules and even functions if possible, I want a map of the project and its distinct areas."
`docs/design/architecture/architecture.html` is the overview, written by hand from `architecture.md`; this page is the code map
under it, and it is generated from the code so it cannot drift from it.

What it reads
=============
Every module of `src/planetgen` (and, as a module level only, `src/tests`) is parsed with the `ast` module. No
module is imported or run. The page is one self-contained file with no network and no build step.

The levels, each reached by clicking down from the one above (every view has a stable `#/...` anchor):

1. **The project map** (`#/`): the areas (web, api, queue, generation, galaxy, population, physics, describe, db,
   admin, backup, cli, names, util, tuning and the tests) as boxes, with the imports between them and the
   purpose line of each package (the first line of its `__init__.py` docstring, which `architecture.md` repeats).
2. **An area** (`#/area/<area>`): its modules, what each imports and what imports it.
3. **A module** (`#/mod/<module>`): its classes and functions with the first line of each docstring, and for each
   function the functions it calls and the functions that call it.
4. **A function** (`#/fn/<module>/<name>`): its source, with every call edge a link.
5. **The flows** (`#/flow`): the main paths end to end (`FLOWS` below, a short hand-kept list of entry
   functions); the generator traces each chain through the resolved calls, so the chain is read from the code
   and nothing in it is drawn by hand.

How a call is resolved
======================
A call is resolved by what the module itself says: a function or class defined in it, a name or module it imports
(`from a import b`, `import a.b as c`, followed through a package's own re-exports), `self.method` and `cls.method`
in the class (and the base classes it names), and `Class.method`. Those are shown as plain links. A method call on
a value of unknown type (`thing.save()`) is resolved only when exactly one function or method in the whole
project has that name (and it is not a name every container has, such as `get`), and is then shown as "probable";
every other call is listed as unresolved, by its text, so the reader sees that it exists and that the map could not
say where it goes. Calls to builtins are left out.

The committed page and the check
================================
The page carries a "Last updated" banner (the Eastern time of the last rebuild that changed it).
The page is committed, like the TODO pages. `--check` fails when a module of `src/planetgen` is missing from the
page or the page names a module that is gone (the list of modules sits in the page as plain JSON), so a lane that
adds or moves a module regenerates the page (test modules are on the page but not checked, so adding a test
file costs no regeneration); on a merge conflict in the page, take main's and run the script
again. The imports, functions and calls are checked only by `--check --full`, which the release stamp can run: the
module list is what a lane's change must keep right, and the rest follows from rerunning the script.

Usage
=====
    python scripts/build_architecture_map.py            # write docs/generated/architecture-map.html
    python scripts/build_architecture_map.py --check    # fail if the module list is out of date
    python scripts/build_architecture_map.py --check --full   # fail if anything differs
"""

import argparse
import ast
import base64
import builtins
import datetime
import json
import os
import re
import sys
import zlib

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src")
PACKAGE = "planetgen"
PAGE = os.path.join(ROOT, "docs", "generated", "architecture-map.html")

AREA_ORDER = ("web", "api", "queue", "generation", "galaxy", "population", "physics", "describe", "db", "admin",
              "backup", "cli", "datasets", "names", "util", "tuning", "package", "tests")
"""tuple[str]: The areas, in the order the project map and the lists show them. An area is the first folder under
`src/planetgen`; `tuning.py` is its own area and so are the package's two root files (`package`); `src/tests` is
`tests`."""

AREA_PURPOSE = {
    "tuning": "Constants specific to this program's design, not physical fact: the generation-tuning knobs and "
              "probabilities.",
    "package": "The package root: the version stamped by the release workflow.",
    "tests": "The test suite (modules and what they import; functions are not mapped).",
}
"""dict: The purpose line of an area that has no package `__init__.py` docstring."""

FLOWS = [
    {"key": "generate", "title": "Generate a galaxy from a command",
     "text": "`planetgen plan` builds the galaxy's skeleton and scatters; `planetgen galaxy` fills sectors, links "
             "them and settles their paths. Each stage runs its work through the work queue and writes through "
             "`db/store`.",
     "entries": ["planetgen.generation.run_plan:run_plan", "planetgen.generation.run_galaxy:run_galaxy"]},
    {"key": "scatter", "title": "A scatter pass (stars and phenomena)",
     "text": "The bright-star scatter and the phenomena scatter: shells of the galaxy are drawn in parallel and "
             "the rows are written to the database.",
     "entries": ["planetgen.generation.run_plan:scatter_bright_stars",
                 "planetgen.generation.run_plan:scatter_phenomena"]},
    {"key": "page", "title": "A page request, from the web route to the database",
     "text": "A browser asks for a page; the route builds its model from the database and renders it.",
     "entries": ["planetgen.web.sector_page:sector"]},
    {"key": "api-job", "title": "An API or Generate-page job through the work queue",
     "text": "A request that starts work is queued and run by a worker, which calls the same generation code the "
             "command line does.",
     "entries": ["planetgen.api.routes:create_system_job", "planetgen.queue.api_jobs:submit"]},
    {"key": "node", "title": "A map node read",
     "text": "A Galaxy Map node request reads the stored node or builds it from the database.",
     "entries": ["planetgen.api.routes:galaxy_nodes_route"]},
    {"key": "backup", "title": "A backup",
     "text": "The command line and the admin page start a database backup, write it and list it.",
     "entries": ["planetgen.cli.backup:run_take", "planetgen.backup.core:take"]},
]
"""list[dict]: The flow view's chains. Each is an entry function (`module:qualified.name`) or several; the
generator traces what they call. A name that is not in the code is an error, so the list cannot rot."""

FLOW_DEPTH = 5
"""int: How many calls deep the flow view traces from an entry before it stops."""

COMMON_METHODS = frozenset("""
    get set pop add append extend insert remove clear copy update items keys values index count sort reverse join
    split rsplit strip lstrip rstrip replace format encode decode startswith endswith find lower upper title read
    write close open flush seek tell readline readlines writelines run start stop send recv put wait cancel
    discard union intersection difference setdefault fetchone fetchall fetchmany execute executemany commit
    rollback cursor connect result done submit map filter reduce next iter exit enter acquire release is_set
    group groups match search sub compile mean sum min max abs round
""".split())
"""frozenset[str]: Method names that most containers, files and library objects have: a call to one of them on an
unknown value is never resolved by name."""

BUILTIN_NAMES = frozenset(dir(builtins))


# -- reading the code ----------------------------------------------------------------------------------------------

def _first_line(doc):
    if not doc:
        return ""
    for line in doc.strip().splitlines():
        line = line.strip()
        if line:
            return line
    return ""


def _module_files():
    """`{dotted name: path}` of every module of the package and of the tests."""
    found = {}
    for base, prefix in ((os.path.join(SRC, PACKAGE), PACKAGE), (os.path.join(SRC, "tests"), "tests")):
        for folder, dirs, files in os.walk(base):
            dirs[:] = sorted(d for d in dirs if d != "__pycache__")
            for name in sorted(files):
                if not name.endswith(".py"):
                    continue
                relative = os.path.relpath(os.path.join(folder, name), base)
                parts = relative[:-3].split(os.sep)
                if parts[-1] == "__init__":
                    parts = parts[:-1]
                dotted = ".".join([prefix] + parts)
                found[dotted] = os.path.join(folder, name)
    return found


def area_of(module):
    parts = module.split(".")
    if parts[0] == "tests":
        return "tests"
    if len(parts) == 1:
        return "package"
    if parts[1] == "tuning":
        return "tuning"
    if parts[1] == "_version":
        return "package"
    return parts[1]


class _Module:
    """One parsed module: its docstring, imports, classes and functions."""

    def __init__(self, name, path, source, is_package):
        self.name, self.path, self.source, self.is_package = name, path, source, is_package
        self.tree = ast.parse(source, filename=path)
        self.doc = _first_line(ast.get_docstring(self.tree))
        self.imports = set()            # dotted modules this one imports (resolved later against the index)
        self.from_names = {}            # local name -> (dotted module, attribute or None)
        self.module_aliases = {}        # local name -> dotted module
        self.defs = {}                  # top-level name -> ("function", qualname) or ("class", qualname)
        self.classes = {}               # class name -> {"methods": {name: qualname}, "bases": [names]}
        self.functions = []             # (qualname, kind, node, class name or None)

    def package_of(self):
        return self.name if self.is_package else self.name.rpartition(".")[0]


def _parse_module(name, path, index_names):
    with open(path, encoding="utf-8") as handle:
        source = handle.read()
    module = _Module(name, path, source, os.path.basename(path) == "__init__.py")
    for node in module.tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            module.defs[node.name] = ("function", node.name)
            module.functions.append((node.name, "function", node, None))
        elif isinstance(node, ast.ClassDef):
            info = {"methods": {}, "bases": [_dotted(base) for base in node.bases if _dotted(base)]}
            module.classes[node.name] = info
            module.defs[node.name] = ("class", node.name)
            module.functions.append((node.name, "class", node, None))
            for member in node.body:
                if isinstance(member, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    qual = f"{node.name}.{member.name}"
                    info["methods"][member.name] = qual
                    module.functions.append((qual, "method", member, node.name))
    for node in ast.walk(module.tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                target = alias.name
                module.imports.add(target)
                if alias.asname:
                    module.module_aliases[alias.asname] = target
                else:
                    module.module_aliases[target.split(".")[0]] = target.split(".")[0]
        elif isinstance(node, ast.ImportFrom):
            base = node.module or ""
            if node.level:
                parts = module.package_of().split(".")
                parts = parts[:len(parts) - (node.level - 1)] if node.level > 1 else parts
                base = ".".join(parts + ([node.module] if node.module else []))
            module.imports.add(base)
            for alias in node.names:
                local = alias.asname or alias.name
                if alias.name == "*":
                    continue
                if f"{base}.{alias.name}" in index_names:
                    module.imports.add(f"{base}.{alias.name}")
                    module.module_aliases[local] = f"{base}.{alias.name}"
                else:
                    module.from_names[local] = (base, alias.name)
    return module


def _dotted(node):
    """`a.b.c` of a Name or Attribute chain, else `None`."""
    parts = []
    while isinstance(node, ast.Attribute):
        parts.append(node.attr)
        node = node.value
    if isinstance(node, ast.Name):
        parts.append(node.id)
        return ".".join(reversed(parts))
    return None


# -- resolving calls -----------------------------------------------------------------------------------------------

class _Resolver:
    def __init__(self, modules):
        self.modules = modules
        self.by_method_name = {}
        for module in modules.values():
            if module.name.startswith("tests"):
                continue
            for qual, kind, _node, _cls in module.functions:
                if kind in ("function", "method"):
                    self.by_method_name.setdefault(qual.rpartition(".")[2], []).append(f"{module.name}:{qual}")

    def lookup(self, module_name, attr, depth=0):
        """The node id of `attr` in `module_name` (a function, or a class's constructor), following the module's own
        re-exports a few steps; `None` when it is not there."""
        module = self.modules.get(module_name)
        if module is None or depth > 3:
            return None
        if attr in module.defs:
            return self._node(module, attr)
        if attr in module.from_names:
            target_module, target_attr = module.from_names[attr]
            return self.lookup(target_module, target_attr, depth + 1)
        return None

    def _node(self, module, name):
        kind, qual = module.defs[name]
        if kind == "class":
            methods = module.classes[name]["methods"]
            if "__init__" in methods:
                return f"{module.name}:{methods['__init__']}"
        return f"{module.name}:{qual}"

    def method_of(self, module, class_name, attr, depth=0):
        """`class_name.attr` of `module`, or of a base class of it."""
        info = module.classes.get(class_name)
        if info is None or depth > 4:
            return None
        if attr in info["methods"]:
            return f"{module.name}:{info['methods'][attr]}"
        for base in info["bases"]:
            head = base.split(".")[-1]
            if head in module.classes:
                found = self.method_of(module, head, attr, depth + 1)
                if found:
                    return found
            elif base.split(".")[0] in module.from_names:
                target_module, target_attr = module.from_names[base.split(".")[0]]
                other = self.modules.get(target_module)
                if other is not None:
                    found = self.method_of(other, target_attr, attr, depth + 1)
                    if found:
                        return found
        return None

    def resolve(self, module, class_name, call):
        """`(kind, value)` of a call: `("exact", node id)`, `("probable", node id)`, `("unresolved", text)` or
        `None` for a builtin."""
        func = call.func
        if isinstance(func, ast.Name):
            name = func.id
            if name in module.defs:
                return "exact", self._node(module, name)
            if name in module.from_names:
                found = self.lookup(*module.from_names[name])
                if found:
                    return "exact", found
            if name in BUILTIN_NAMES:
                return None
            return "unresolved", name
        parts = _dotted(func)
        if parts is None:
            if isinstance(func, ast.Attribute):
                return self._by_name(func.attr, "." + func.attr)
            return "unresolved", "(...)"
        pieces = parts.split(".")
        head = pieces[0]
        if head in ("self", "cls") and class_name and len(pieces) == 2:
            found = self.method_of(module, class_name, pieces[1])
            if found:
                return "exact", found
            return self._by_name(pieces[1], parts)
        if head in module.module_aliases:
            target = module.module_aliases[head]
            rest = pieces[1:]
            while len(rest) > 1 and f"{target}.{rest[0]}" in self.modules:
                target, rest = f"{target}.{rest[0]}", rest[1:]
            if len(rest) == 1:
                found = self.lookup(target, rest[0])
                if found:
                    return "exact", found
            elif len(rest) == 2 and target in self.modules:
                other = self.modules[target]
                if rest[0] in other.classes:
                    found = self.method_of(other, rest[0], rest[1])
                    if found:
                        return "exact", found
            return "unresolved", parts
        if head in module.classes and len(pieces) == 2:
            found = self.method_of(module, head, pieces[1])
            if found:
                return "exact", found
        if head in module.from_names and len(pieces) == 2:
            target_module, target_attr = module.from_names[head]
            other = self.modules.get(target_module)
            if other is not None and target_attr in other.classes:
                found = self.method_of(other, target_attr, pieces[1])
                if found:
                    return "exact", found
        return self._by_name(pieces[-1], parts)

    def _by_name(self, attr, text):
        if attr not in COMMON_METHODS and not attr.startswith("__"):
            found = self.by_method_name.get(attr, [])
            if len(found) == 1:
                return "probable", found[0]
        return "unresolved", text


def _calls_in(node):
    for child in ast.walk(node):
        if isinstance(child, ast.Call):
            yield child


# -- the data ------------------------------------------------------------------------------------------------------

def collect():
    """
    Reads the whole code base.

    Returns:
        dict: `modules` (name -> `_Module`), `files` (name -> path).
    """
    files = _module_files()
    names = set(files)
    modules = {name: _parse_module(name, path, names) for name, path in files.items()}
    return {"modules": modules, "files": files}


def build_data(collected=None):
    """
    The page's data: areas, modules with their imports, functions with their calls, the flows and the sources.

    Returns:
        tuple[dict, dict]: `(light, heavy)`: `light` is the module list and areas (written in the page as plain
        JSON, what `--check` compares); `heavy` is the rest (written compressed).
    """
    collected = collected or collect()
    modules = collected["modules"]
    resolver = _Resolver(modules)
    order = {area: n for n, area in enumerate(AREA_ORDER)}
    module_names = sorted(modules, key=lambda m: (order.get(area_of(m), 99), m))
    module_index = {name: n for n, name in enumerate(module_names)}

    areas = []
    for area in AREA_ORDER:
        members = [m for m in module_names if area_of(m) == area]
        if not members:
            continue
        init = modules.get(f"{PACKAGE}.{area}")
        purpose = AREA_PURPOSE.get(area) or (init.doc if init is not None else "")
        areas.append({"key": area, "purpose": purpose, "modules": [module_index[m] for m in members]})

    light = {"areas": areas,
             "modules": [{"name": name, "area": area_of(name), "doc": modules[name].doc,
                          "package": modules[name].is_package} for name in module_names]}

    functions, function_index = [], {}
    for name in module_names:
        if name.startswith("tests"):
            continue
        for qual, kind, node, class_name in modules[name].functions:
            function_index[f"{name}:{qual}"] = len(functions)
            functions.append({"id": f"{name}:{qual}", "module": module_index[name], "name": qual, "kind": kind,
                              "doc": _first_line(ast.get_docstring(node)), "line": node.lineno,
                              "end": node.end_lineno, "node": node, "class": class_name})

    imports = []
    for name in module_names:
        found = sorted({module_index[target] for target in modules[name].imports if target in modules
                        and target != name})
        imports.append(found)

    packed_functions, unresolved_total = [], 0
    for entry in functions:
        module = modules[module_names[entry["module"]]]
        calls, text = [], []
        seen = set()
        if entry["kind"] != "class":
            for call in _calls_in(entry["node"]):
                result = resolver.resolve(module, entry["class"], call)
                if result is None:
                    continue
                kind, value = result
                if kind == "unresolved":
                    if value not in seen:
                        seen.add(value)
                        text.append(value)
                    continue
                target = function_index.get(value)
                if target is None or value == entry["id"] or (target, kind) in seen:
                    continue
                seen.add((target, kind))
                calls.append([target, 0 if kind == "exact" else 1])
        unresolved_total += len(text)
        packed_functions.append([entry["module"], entry["name"], entry["kind"], entry["doc"], entry["line"],
                                 entry["end"], calls, text])

    flows = []
    for flow in FLOWS:
        entries = []
        for entry in flow["entries"]:
            if entry not in function_index:
                raise SystemExit(f"flow {flow['key']!r}: no function {entry!r} in the code; fix FLOWS in "
                                 f"scripts/build_architecture_map.py")
            entries.append(function_index[entry])
        flows.append({"key": flow["key"], "title": flow["title"], "text": flow["text"], "entries": entries})

    sources = {module_index[name]: modules[name].source for name in module_names if not name.startswith("tests")}
    heavy = {"imports": imports, "functions": packed_functions, "flows": flows, "depth": FLOW_DEPTH,
             "sources": {str(k): v for k, v in sources.items()}}
    return light, heavy


# -- the page ------------------------------------------------------------------------------------------------------

TEMPLATE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "architecture_map", "template.html")


def render(light, heavy):
    """The page's text."""
    with open(TEMPLATE_PATH, encoding="utf-8") as handle:
        template = handle.read()
    blob = base64.b64encode(zlib.compress(json.dumps(heavy, separators=(",", ":")).encode("utf-8"), 9)).decode("ascii")
    light_text = json.dumps(light, separators=(",", ":")).replace("</", "<\\/")
    return template.replace("/*LIGHT*/", light_text).replace("/*HEAVY*/", blob)


STAMP_RE = re.compile(r'<div class="updated-banner"[^>]*>Last updated: (.*?)</div>\n')


def with_stamp(page, stamp):
    """`page` with a "Last updated" banner (Boss, 3:50 PM EDT 2026-10-31) just above its header."""
    style = ("font:600 13px/1.4 system-ui,-apple-system,'Segoe UI',sans-serif;padding:6px 14px;text-align:center;"
             "background:#F6E7D6;color:#152030;border-bottom:1px solid #D8DEE7")
    page = STAMP_RE.sub("", page)
    banner = f'<div class="updated-banner" style="{style}">Last updated: {stamp}</div>\n'
    at = page.index("<header>")
    return page[:at] + banner + page[at:]


def now_stamp():
    """The Eastern time now, for example `2026-10-31 4:05 PM EDT`."""
    from zoneinfo import ZoneInfo
    moment = datetime.datetime.now(ZoneInfo("America/New_York"))
    return f"{moment:%Y-%m-%d} {moment.hour % 12 or 12}:{moment:%M %p} {moment:%Z}"


def build(stamp=None):
    """The page. The banner keeps the stamp the committed page carries (so `--check --full` stays exact); `main`
    writes a new time only when the rest of the page changed."""
    if stamp is None:
        found = STAMP_RE.search(open(PAGE, encoding="utf-8").read()) if os.path.exists(PAGE) else None
        stamp = found.group(1) if found else ""
    return with_stamp(render(*build_data()), stamp)


LIGHT_RE = re.compile(r'<script type="application/json" id="light">(.*?)</script>', re.S)


def page_light(text):
    """The module list a page text carries, `None` if it has none."""
    found = LIGHT_RE.search(text)
    return None if found is None else json.loads(found.group(1).replace("<\\/", "</"))


def module_problems(text, collected=None):
    """What is wrong with the page's module list against the code: modules missing from the page and modules the
    page names that are gone (`[]` when it is right)."""
    light = page_light(text)
    if light is None:
        return ["the page has no module list; run python scripts/build_architecture_map.py"]
    on_page = {m["name"] for m in light["modules"] if m["area"] != "tests"}
    in_code = {name for name in (collected or collect())["files"] if area_of(name) != "tests"}
    problems = [f"missing from the page: {name}" for name in sorted(in_code - on_page)]
    problems += [f"on the page but not in the code: {name}" for name in sorted(on_page - in_code)]
    return problems


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--check", action="store_true", help="fail if the page's module list is out of date")
    parser.add_argument("--full", action="store_true", help="with --check: fail if anything on the page differs")
    args = parser.parse_args(argv)
    if args.check:
        if not os.path.exists(PAGE):
            print("docs/generated/architecture-map.html is missing; run python scripts/build_architecture_map.py")
            return 1
        with open(PAGE, encoding="utf-8") as handle:
            current = handle.read()
        problems = module_problems(current)
        if args.full and not problems and current != build():
            problems = ["the page differs from the code; run python scripts/build_architecture_map.py"]
        for problem in problems:
            print(problem)
        if problems:
            print("docs/generated/architecture-map.html is out of date: run python scripts/build_architecture_map.py")
            return 1
        print("architecture map is up to date")
        return 0
    text = build()
    old = open(PAGE, encoding="utf-8").read() if os.path.exists(PAGE) else None
    if text != old:
        text = build(now_stamp())
    with open(PAGE, "w", encoding="utf-8") as handle:
        handle.write(text)
    light = page_light(text)
    print(f"wrote {os.path.relpath(PAGE, ROOT)}: {len(light['modules'])} modules in {len(light['areas'])} areas, "
          f"{len(text) / 1e6:.2f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
