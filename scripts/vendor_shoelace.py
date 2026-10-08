"""Copy the parts of Shoelace the site uses into src/html/static/vendor/shoelace/.

Shoelace (https://shoelace.style, MIT) ships its `cdn/` folder as native ES
modules: one file per component plus shared `chunks/`. The pages load them
straight from the site's own `static/` (the Content-Security-Policy is
`default-src 'self'`, and there is no bundler). This script follows each
wanted component's `import` lines and copies only the files they reach, so
the vendored set stays small.

Usage:

    npm pack @shoelace-style/shoelace@<version>
    tar xzf shoelace-style-shoelace-<version>.tgz
    python scripts/vendor_shoelace.py package

To use another component, add it to COMPONENTS and run the script again;
add icons to ICONS the same way. It rewrites the whole destination folder.
"""

import json
import os
import re
import shutil
import subprocess
import sys
import urllib.parse

COMPONENTS = ("button", "icon-button", "icon", "dropdown", "menu", "menu-item", "divider")
"""tuple: Shoelace components the shared set (static/components.js) registers."""

UTILITIES = ("icon-library",)
"""tuple: Shoelace utility modules the shared set imports (static/shoelaceicons.js)."""

ICONS = ("gear",)
"""tuple: Bootstrap icons (Shoelace's icon set, MIT) copied to assets/icons/."""

SYSTEM_ICONS = ("caret", "check", "chevron-left", "chevron-right")
"""tuple: Icons Shoelace's components draw themselves (the caret on a button,
the check in a menu item), written to assets/system/. Shoelace holds them as
`data:` URIs, which the site's Content-Security-Policy refuses to fetch, so
they are saved as files and static/shoelaceicons.js points the "system"
library at them. Add a name here when a component you start using shows one
(its icon request 404s under assets/system/)."""

DEST = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src", "html", "static", "vendor", "shoelace")

_IMPORT = re.compile(r'(?:from|import)\s*"(\.[^"]+)"|import\("(\.[^"]+)"\)')


def closure(cdn, entries):
    """Return the set of files (relative to `cdn`) reachable from `entries` by relative imports."""
    seen = set()
    todo = [os.path.normpath(e) for e in entries]
    while todo:
        path = todo.pop()
        if path in seen:
            continue
        seen.add(path)
        with open(os.path.join(cdn, path), encoding="utf-8") as handle:
            source = handle.read()
        for match in _IMPORT.finditer(source):
            target = match.group(1) or match.group(2)
            todo.append(os.path.normpath(os.path.join(os.path.dirname(path), target)))
    return seen


def write_system_icons(cdn, dest):
    """Save the SYSTEM_ICONS of `chunks/*` (found by their `name: "system"`) as SVG files; needs node."""
    chunks = os.path.join(cdn, "chunks")
    holder = next(name for name in sorted(os.listdir(chunks))
                  if 'name: "system"' in open(os.path.join(chunks, name), encoding="utf-8").read())
    script = ("import(process.argv[1]).then(m => { const lib = m.library_system_default; "
              "console.log(JSON.stringify(process.argv.slice(2).map(n => lib.resolver(n)))); });")
    out = subprocess.run(["node", "-e", script, os.path.abspath(os.path.join(chunks, holder)), *SYSTEM_ICONS],
                         check=True, capture_output=True, text=True).stdout
    os.makedirs(os.path.join(dest, "assets", "system"), exist_ok=True)
    for name, uri in zip(SYSTEM_ICONS, json.loads(out)):
        svg = urllib.parse.unquote(uri.split(",", 1)[1])
        with open(os.path.join(dest, "assets", "system", f"{name}.svg"), "w", encoding="utf-8") as handle:
            handle.write(svg.strip() + "\n")


def main(package):
    cdn = os.path.join(package, "cdn")
    dest = os.path.normpath(DEST)
    shutil.rmtree(dest, ignore_errors=True)
    entries = [f"components/{name}/{name}.js" for name in COMPONENTS]
    entries += [f"utilities/{name}.js" for name in UTILITIES]
    files = closure(cdn, entries)
    files.add("themes/light.css")
    for icon in ICONS:
        files.add(f"assets/icons/{icon}.svg")
    for path in sorted(files):
        target = os.path.join(dest, path)
        os.makedirs(os.path.dirname(target), exist_ok=True)
        shutil.copyfile(os.path.join(cdn, path), target)
    write_system_icons(cdn, dest)
    shutil.copyfile(os.path.join(package, "LICENSE.md"), os.path.join(dest, "LICENSE.md"))
    print(f"{len(files)} files copied to {dest}")


if __name__ == "__main__":
    main(sys.argv[1])
