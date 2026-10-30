#!/usr/bin/env python3
"""
The wiring guard (GEN.754): what the app builds and downloads but nothing reads.

    python scripts/wiring_guard.py            # print every unread thing and whether the register names it
    python scripts/wiring_guard.py --check    # exit 1 when the register and the code disagree

Four kinds of thing can be built and left unread, and the guard finds each from the source alone (no import of the app, so it
runs without the datasets or a database):

- a `module`: a module of `physics`, `generation`, `galaxy`, `datasets` or `tuning` that no module outside those packages
  reaches through imports (the app is the CLI, API, web, admin, queue, population, backup and describe packages; lazy imports
  inside functions count);
- a `dataset`: a dataset of the setup catalogue (`planetgen.datasets.catalogue.ALL`) whose module nothing reachable outside
  `datasets` imports, or a `reader`, a function of such a module that opens the downloaded files and that no reachable module
  calls;
- a `data` file: a `physics/data/*.json` whose file name no reachable module mentions;
- a `row`: a registry row (GEN.725, `gap-filling-plan.md` section 6) whose `module.CONSTANT` names are all unread, so a stand-in
  that nothing in the generator uses.

`docs/design/research/wiring-register.md` names every unread thing with the open item that will wire it (or deletes it), or says
why it is read by design (a tuning harness loaded by name). `src/tests/test_wiring_guard.py` fails when something is unread and
not in the register, when a register line names something that is read now, and when a named item is not open in
`docs/TODO.md`. No shim: an unread piece is wired, filed or deleted.
"""

import argparse
import ast
import collections
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
SRC = os.path.join(ROOT, "src")
PACKAGE = os.path.join(SRC, "planetgen")
REGISTER = os.path.join(ROOT, "docs", "design", "research", "wiring-register.md")
TODO = os.path.join(ROOT, "docs", "TODO.md")
PLAN = os.path.join(ROOT, "docs", "design", "research", "gap-filling-plan.md")
DATA_FOLDER = os.path.join(PACKAGE, "physics", "data")

LIBRARY = ("planetgen.physics", "planetgen.generation", "planetgen.galaxy", "planetgen.datasets", "planetgen.tuning")
"""tuple: The packages and modules whose pieces must be reached from the rest of the app."""

KINDS = ("module", "dataset", "reader", "data", "row")
BY_DESIGN = "by design"
"""str: The item cell of a thing that is read by name or only by the tuning harness: the reason says which."""

_CONSTANT = re.compile(r"\b([a-z_][a-z_0-9]*)\.([A-Z][A-Z0-9_]{2,})\b")


def _tree(path):
    with open(path, encoding="utf-8") as handle:
        return ast.parse(handle.read(), filename=path)


def modules():
    """Every module of the package, as `{dotted name: path}`."""
    found = {}
    for folder, names, files in os.walk(PACKAGE):
        names[:] = [name for name in names if name != "__pycache__"]
        for name in files:
            if name.endswith(".py"):
                path = os.path.join(folder, name)
                dotted = os.path.relpath(path, SRC)[:-3].replace(os.sep, ".")
                found[dotted[:-9] if dotted.endswith(".__init__") else dotted] = path
    return found


def _imports(name, path, known, tree):
    """The modules of the package that `tree` imports, with the parents Python imports first."""
    package = name if path.endswith("__init__.py") else name.rpartition(".")[0]
    wanted = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            wanted.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                parts = package.split(".")
                module = ".".join(parts[:len(parts) - node.level + 1] + ([node.module] if node.module else []))
            else:
                module = node.module or ""
            wanted.add(module)
            wanted.update(f"{module}.{alias.name}" for alias in node.names)
    found = set()
    for dotted in wanted:
        parts = dotted.split(".")
        found.update(".".join(parts[:i]) for i in range(1, len(parts) + 1) if ".".join(parts[:i]) in known)
    return found


class Code:
    """The parsed package: modules, who imports whom, what is reachable from the app, and the names it uses."""

    def __init__(self):
        self.paths = modules()
        self.trees = {name: _tree(path) for name, path in self.paths.items()}
        self.imports = {name: _imports(name, self.paths[name], self.paths, tree) for name, tree in self.trees.items()}
        self.named = {name: self._named_modules(tree) for name, tree in self.trees.items()}
        self.reachable = self._reachable()
        self.loaded = collections.defaultdict(set)
        self.strings = collections.defaultdict(set)
        for name in self.reachable:
            for node in ast.walk(self.trees[name]):
                if isinstance(node, ast.Attribute):
                    self.loaded[node.attr].add(name)
                elif isinstance(node, ast.Name) and isinstance(node.ctx, ast.Load):
                    self.loaded[node.id].add(name)
                elif isinstance(node, ast.ImportFrom):
                    for alias in node.names:
                        self.loaded[alias.name].add(name)
                elif isinstance(node, ast.Constant) and isinstance(node.value, str):
                    self.strings[node.value].add(name)

    def _named_modules(self, tree):
        """The modules of the package that a string of `tree` names (`importlib.import_module("planetgen.x")`, or the
        `"planetgen.x:function"` form of a registry that loads its builders by name)."""
        found = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and isinstance(node.value, str) and node.value.startswith("planetgen."):
                dotted = node.value.split(":")[0]
                if dotted in self.paths:
                    found.add(dotted)
        return found

    def _reachable(self):
        roots = [name for name in self.paths if not name.startswith(LIBRARY)]
        seen, stack = set(roots), list(roots)
        while stack:
            name = stack.pop()
            for other in self.imports[name] | self.named[name]:
                if other not in seen:
                    seen.add(other)
                    stack.append(other)
        return seen

    def unread_modules(self):
        return sorted(name for name in self.paths if name.startswith(LIBRARY) and name not in self.reachable)

    def _outside_importers(self, module):
        return sorted(name for name in self.reachable if module in self.imports[name] and not name.startswith("planetgen.datasets"))

    def unread_datasets(self):
        """`(dataset module names that nothing reachable imports, 'module.function' readers that nothing reachable calls)`."""
        catalogue = self.trees["planetgen.datasets.catalogue"]
        names = []
        for node in ast.walk(catalogue):
            if isinstance(node, ast.ImportFrom) and node.module == "planetgen.datasets":
                names += [alias.name for alias in node.names if alias.name != "fetch"]
        datasets, readers = [], []
        for short in sorted(names):
            module = f"planetgen.datasets.{short}"
            users = self._outside_importers(module)
            if not users:
                datasets.append(short)
                continue
            used = set()
            for user in users:
                for node in ast.walk(self.trees[user]):
                    if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name) and node.value.id == short:
                        used.add(node.attr)
                    elif isinstance(node, ast.ImportFrom) and node.module == module:
                        used.update(alias.name for alias in node.names)
            for node in self.trees[module].body:
                if isinstance(node, ast.FunctionDef) and not node.name.startswith("_") and node.name not in used \
                        and self._opens_files(node):
                    readers.append(f"{short}.{node.name}")
        return datasets, readers

    @staticmethod
    def _opens_files(function):
        return any(isinstance(node, ast.Name) and node.id in ("fetch", "DATASET") for node in ast.walk(function))

    def unread_data_files(self):
        found = []
        for name in sorted(os.listdir(DATA_FOLDER)):
            if name.endswith(".json") and not any(name in text or text == name[:-5] for text in self.strings):
                found.append(name)
        return found

    def unread_rows(self, tables):
        """The registry rows `{section}:{fill_id}` whose `module.CONSTANT` names all lack a reader in the reachable code."""
        by_base = collections.defaultdict(list)
        for name in self.paths:
            by_base[name.rpartition(".")[2]].append(name)
        found = []
        for section, (header, rows) in tables.items():
            lowered = [cell.lower() for cell in header]
            if "constant" not in lowered:
                continue
            for row in rows:
                named = [constant for module, constant in _CONSTANT.findall(row[lowered.index("constant")])
                         if len(by_base.get(module, ())) == 1]
                if named and not any(self.loaded.get(constant) for constant in named):
                    found.append(f"{section}:{row[0]}")
        return found


def unread(code=None, tables=None):
    """Every unread thing as `{(kind, name)}`."""
    sys.path.insert(0, HERE)
    import build_stand_in_values as builder  # noqa: E402
    code = code or Code()
    if tables is None:
        with open(PLAN, encoding="utf-8") as handle:
            tables = builder.registry_tables(handle.read())
    datasets, readers = code.unread_datasets()
    found = {("module", name) for name in code.unread_modules()}
    found |= {("dataset", name) for name in datasets} | {("reader", name) for name in readers}
    found |= {("data", name) for name in code.unread_data_files()}
    found |= {("row", name) for name in code.unread_rows(tables)}
    return found


def parse_register(text):
    """The register's lines as `{(kind, name): (item, reason)}`; `item` is `by design` or an item id such as `GEN.850`."""
    entries = {}
    for line in text.splitlines():
        cells = [cell.strip().strip("`") for cell in line.strip().strip("|").split("|")]
        if line.startswith("|") and len(cells) == 4 and cells[0] in KINDS:
            entries[(cells[0], cells[1])] = (cells[2], cells[3])
    return entries


def load_register():
    with open(REGISTER, encoding="utf-8") as handle:
        return parse_register(handle.read())


def open_items(text):
    """The ids of the open items of `docs/TODO.md`."""
    return set(re.findall(r"^- \[ \] \*\*([A-Z]+\.\d+[a-z]?)\b", text, re.MULTILINE))


def problems(found, entries, todo_text):
    """The disagreements between the code (`found`) and the register (`entries`), one sentence each."""
    opened = open_items(todo_text)
    messages = []
    for key in sorted(found - set(entries)):
        messages.append(f"{key[0]} {key[1]} is built and nothing reads it: wire it, file an open item for it and name it in the "
                        f"register, or delete it")
    for key in sorted(set(entries) - found):
        messages.append(f"the register names {key[0]} {key[1]} but it is read now (or gone): delete the line")
    for key, (item, reason) in sorted(entries.items()):
        if not reason:
            messages.append(f"{key[0]} {key[1]}: the register line gives no reason")
        if item != BY_DESIGN and not (set(re.findall(r"[A-Z]+\.\d+[a-z]?", item)) & opened):
            messages.append(f"{key[0]} {key[1]}: item {item!r} is not an open item of docs/TODO.md")
    return messages


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("--check", action="store_true", help="exit 1 when the register and the code disagree")
    args = parser.parse_args(argv)
    found, entries = unread(), load_register()
    with open(TODO, encoding="utf-8") as handle:
        messages = problems(found, entries, handle.read())
    for key in sorted(found):
        item = entries.get(key, ("NOT IN THE REGISTER", ""))[0]
        print(f"{key[0]:8} {key[1]:60} {item}")
    for message in messages:
        print("PROBLEM:", message, file=sys.stderr)
    return 1 if args.check and messages else 0


if __name__ == "__main__":
    sys.exit(main())
