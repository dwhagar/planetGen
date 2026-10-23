# tests/test_unit_registry_js.py

"""
UX.140: the unit registry (`planetgen.physics.unit_registry`) is written to the browser as the generated
`html/static/unitregistry.js`. The file is current, every unit's numbers and every conversion in it give exactly what the
Python registry gives (the same arithmetic on the same doubles, checked through Node), and no browser script declares a
unit factor of its own.
"""

import itertools
import json
import os
import re
import shutil
import subprocess
import sys

import pytest

from planetgen.physics import unit_registry as registry

_ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
_STATIC = os.path.join(_ROOT, "src", "html", "static")
_MODULE = os.path.join(_STATIC, "unitregistry.js")

_SAMPLES = (0.0, 1.0, -1.0, 1.5, 0.001, 3.7, 12345.678, 6.02e23, 1.2e-9, 299792458.0, 273.15, -40.0, 1e-30, 1e30)

_NODE = pytest.mark.skipif(shutil.which("node") is None, reason="Node isn't installed")

_SCRIPT = (
    "import fs from 'node:fs';"
    "const m = await import(process.argv[1]);"
    "const job = JSON.parse(fs.readFileSync(0, 'utf8'));"
    "const out = {"
    "  units: Object.fromEntries(Object.entries(m.UNITS).map(([k, u]) => [k, [u.kind, u.symbol, u.name, u.singular, u.factor, u.offset]])),"
    "  kinds: Object.fromEntries(Object.entries(m.KINDS).map(([k, v]) => [k, [v.base, v.dimension, v.low, v.high, v.lowOpen, v.highOpen]])),"
    "  converted: job.conversions.map(([v, a, b]) => m.convert(v, a, b)),"
    "  based: job.bases.map(([v, k]) => [m.toBase(v, k), m.fromBase(v, k)]),"
    "};"
    "console.log(JSON.stringify(out));")


def _pairs():
    for kind in registry.KINDS:
        for source, target in itertools.permutations([unit.key for unit in registry.units_of(kind)], 2):
            yield kind, source, target


def _in_node(job):
    run = subprocess.run(["node", "--input-type=module", "-e", _SCRIPT, "file://" + os.path.abspath(_MODULE)],
                         input=json.dumps(job), capture_output=True, text=True, check=True, timeout=120)
    return json.loads(run.stdout, parse_int=float)


def test_the_committed_file_is_current():
    run = subprocess.run([sys.executable, os.path.join(_ROOT, "scripts", "build_unit_registry.py"), "--check"],
                         capture_output=True, text=True)
    assert run.returncode == 0, run.stdout + run.stderr


@_NODE
def test_every_unit_and_every_conversion_is_the_python_registrys():
    pairs = list(_pairs())
    conversions = [[value, source, target] for _kind, source, target in pairs for value in _SAMPLES]
    bases = [[value, key] for key in registry.UNITS for value in _SAMPLES]
    answer = _in_node({"conversions": conversions, "bases": bases})
    assert set(answer["units"]) == set(registry.UNITS) and set(answer["kinds"]) == set(registry.KINDS)
    for key, row in registry.UNITS.items():
        assert answer["units"][key] == [row.kind, row.symbol, row.name, row.singular or row.name, row.factor, row.offset], key
    for name, row in registry.KINDS.items():
        assert answer["kinds"][name] == [row.base, list(row.dimension), row.bounds.low, row.bounds.high,
                                         row.bounds.low_open, row.bounds.high_open], name
    for (value, source, target), got in zip(conversions, answer["converted"]):
        assert got == registry.convert(value, source, target), (value, source, target)
    for (value, key), got in zip(bases, answer["based"]):
        assert got == [registry.to_base(value, key), registry.from_base(value, key)], (value, key)


@_NODE
def test_a_conversion_between_kinds_is_refused_in_the_browser_as_it_is_here():
    script = ("const m = await import(process.argv[1]); try { m.convert(1, 'km', 'kg'); console.log('no'); }"
              " catch (e) { console.log(e.message); }")
    run = subprocess.run(["node", "--input-type=module", "-e", script, "file://" + os.path.abspath(_MODULE)],
                         capture_output=True, text=True, check=True, timeout=60)
    assert "cannot convert km (length) to kg (mass)" in run.stdout


_FACTOR_NAMES = re.compile(r"\b(?:const|let|var|export const)\s+(?:AU_KM|KM_PER_PC|KM_PER_AU|AU_M|KM_M|PARSEC_M|LIGHTYEAR_M"
                           r"|SECONDS_PER_YEAR)\s*=\s*[-+0-9.]")
_FACTOR_VALUES = re.compile(r"149597870\.?7?|3\.0856775814913\d*e1[36]|9460730472580800|206264\.8|31557600\b|3\.15576e7")
_GENERATED = ("unitregistry.js", "unitladders.js")


def test_no_browser_script_declares_a_unit_factor():
    found = []
    for name in sorted(os.listdir(_STATIC)):
        if not name.endswith(".js") or name in _GENERATED:
            continue
        with open(os.path.join(_STATIC, name), encoding="utf-8") as handle:
            for number, line in enumerate(handle, 1):
                code = line.split("//")[0]
                if _FACTOR_NAMES.search(code) or _FACTOR_VALUES.search(code):
                    found.append(f"{name}:{number}: {line.strip()}")
    assert not found, "unit factors declared in the browser (read them from unitregistry.js):\n" + "\n".join(found)
