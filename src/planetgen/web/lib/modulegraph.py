# planetgen/web/lib/modulegraph.py

"""
The modules a page's script pulls in, for `<link rel="modulepreload">` hints
(MAP.161).

The map scripts import their siblings with `await import(`./x.js${VERSION_QUERY}`)`
(see `fmt.static_url`: every module has one URL per page), and a module only
starts loading once the one that names it has arrived and run. The Galaxy Map's
tree is some 30 files deep, so the first star waited for the files one round
trip after another (7.9 to 8.1 s on a slow 4G link). A page that names the whole
tree up front lets the browser fetch it in parallel; the later `import()` finds
the module already in the browser's module map.

Only the literal form above is followed; a computed name (the Shoelace
components, the terminal add-ons) is not a hint, and loads when asked for.
"""

import os
import re

from planetgen.web.lib import fmt

_IMPORT_RE = re.compile(r"""import\(\s*`\./([A-Za-z0-9_./-]+\.m?js)\$\{VERSION_QUERY\}`\s*\)""")

_cache = {}
"""dict: `(entry, fingerprint) -> list[str]`; a changed static folder has a new fingerprint."""


def _imports_of(name):
    """The module names `name` (a path under `static/`) imports, in source order, once each."""
    try:
        with open(os.path.join(fmt._STATIC_DIR, name), encoding="utf-8") as handle:
            source = handle.read()
    except OSError:
        return []
    return list(dict.fromkeys(match.group(1) for match in _IMPORT_RE.finditer(source)))


def module_graph(entry):
    """
    Every module under `static/` that `entry` imports, directly or not, nearest first
    (breadth first), not including `entry` itself.

    Args:
        entry (str): A module's path relative to `static/`, e.g. `"galaxymap3d.js"`.

    Returns:
        list[str]: Paths relative to `static/`, ready for `static_url`.
    """
    key = (entry, fmt._static_fingerprint())
    found = _cache.get(key)
    if found is None:
        found = []
        seen = {entry}
        queue = [entry]
        while queue:
            name = queue.pop(0)
            for child in _imports_of(name):
                if child not in seen:
                    seen.add(child)
                    found.append(child)
                    queue.append(child)
        _cache.clear()
        _cache[key] = found
    return list(found)
