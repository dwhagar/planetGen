# tests/test_modulegraph.py

"""
`planetgen/web/lib/modulegraph.py` (MAP.161): the modules a page's script imports,
named in `<link rel="modulepreload">` hints so they load in parallel.
"""

import os
import re

from planetgen.web.lib import fmt, modulegraph


def test_the_galaxy_maps_whole_tree_is_found():
    graph = modulegraph.module_graph("galaxymap3d.js")
    assert "vendor/three.module.min.js" in graph
    assert "galaxystageview.js" in graph and "galaxysector.js" in graph, "a module two imports down"
    assert "systemview3d.js" in graph, "a module three or more imports down"
    assert "galaxymap3d.js" not in graph and len(graph) == len(set(graph))


def test_every_listed_file_exists_and_every_import_of_one_is_listed():
    graph = modulegraph.module_graph("galaxymap3d.js")
    for name in graph:
        assert os.path.isfile(os.path.join(fmt._STATIC_DIR, name)), name
    listed = set(graph) | {"galaxymap3d.js"}
    for name in listed:
        if name.startswith("vendor/"):
            continue
        with open(os.path.join(fmt._STATIC_DIR, name), encoding="utf-8") as handle:
            source = handle.read()
        for match in re.finditer(r"import\(`\./([A-Za-z0-9_./-]+\.m?js)\$\{VERSION_QUERY\}`\)", source):
            assert match.group(1) in listed, f"{name} imports {match.group(1)}, which is not in the preload list"


def test_a_missing_entry_has_no_modules():
    assert modulegraph.module_graph("no-such-module.js") == []


def test_the_graph_is_a_copy_the_caller_can_change():
    first = modulegraph.module_graph("galaxymap3d.js")
    first.clear()
    assert modulegraph.module_graph("galaxymap3d.js")
