#!/usr/bin/env python3
"""
Builds the proton range and stopping power file (GEN.453) from the NIST PSTAR tables.

    python scripts/build_proton_range.py                  # fetch from NIST, write src/planetgen/physics/data/proton_range.json
    python scripts/build_proton_range.py --raw DIR        # read (and save, when missing) the raw PSTAR pages in DIR
    python scripts/build_proton_range.py --check --raw DIR   # build in memory and fail if the committed file differs

The file is `src/planetgen/physics/data/proton_range.json`. PSTAR (NIST, ICRU 49) gives each material's total stopping
power and CSDA range for protons at 133 energies from 1 keV to 10 GeV. The file keeps the CSDA range and the total
stopping power at the 26 energies of `ENERGIES_MEV` (1 MeV to 10 GeV), for the materials of `MATERIALS` (the elements
and compounds the atmosphere and ground chains name). Nothing is fitted: a value is a PSTAR row as printed. The raw pages
are never in the repository (DATA-LICENSE.md); the thinned table is.
"""

import argparse
import json
import os
import re
import sys
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "src", "planetgen", "physics", "data", "proton_range.json")
URL = "https://physics.nist.gov/cgi-bin/Star/ap_table-t.pl"
RETRIEVED = "2026-10-23"

ENERGIES_MEV = (1, 2, 3, 5, 7, 10, 15, 20, 30, 40, 50, 70, 100, 150, 200, 300, 400, 500, 700, 1000, 1500, 2000, 3000,
                5000, 7000, 10000)
"""tuple: The proton kinetic energies, MeV, kept from PSTAR's grid (every one is a row of the PSTAR default grid)."""

# key -> (PSTAR material number, PSTAR's name for it)
MATERIALS = {
    "H": ("001", "Hydrogen"), "He": ("002", "Helium"), "C": ("906", "Graphite (density 1.7 g/cm3)"),
    "N": ("007", "Nitrogen"), "O": ("008", "Oxygen"), "Ne": ("010", "Neon"), "Si": ("014", "Silicon"),
    "Ar": ("018", "Argon"), "Fe": ("026", "Iron"), "Kr": ("036", "Krypton"), "Xe": ("054", "Xenon"),
    "air": ("104", "Air, Dry (near sea level)"), "CO2": ("134", "Carbon Dioxide"), "CH4": ("197", "Methane"),
    "SiO2": ("245", "Silicon Dioxide"), "H2O liquid": ("276", "Water, Liquid"), "H2O vapor": ("277", "Water Vapor"),
}

ROW = re.compile(r"^\s*(\d\.\d{3}E[+-]\d\d)((?:\s+\d\.\d{3,4}E?[+-]?\d*){6})\s*$")


def _fetch(number, raw):
    path = os.path.join(raw, "pstar-%s.html" % number) if raw else None
    if path and os.path.exists(path):
        with open(path, encoding="utf-8", errors="replace") as handle:
            return handle.read()
    data = urllib.parse.urlencode({"prog": "PSTAR", "matno": number, "ShowDefault": "on", "Energies": ""}).encode()
    with urllib.request.urlopen(urllib.request.Request(URL, data=data), timeout=60) as response:
        text = response.read().decode("utf-8", errors="replace")
    if path:
        os.makedirs(raw, exist_ok=True)
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(text)
    return text


def _rows(text):
    """The PSTAR rows of a page as `{energy_mev: (stop_e, stop_n, stop_t, range_c, range_p, detour)}`."""
    rows = {}
    for line in re.sub(r"<br\s*/?>", "\n", text).splitlines():
        match = ROW.match(line)
        if match:
            numbers = [float(match.group(1))] + [float(v) for v in match.group(2).split()]
            rows[numbers[0]] = tuple(numbers[1:])
    return rows


def build(raw=None):
    materials = {}
    for key, (number, name) in MATERIALS.items():
        rows = _rows(_fetch(number, raw))
        picked = []
        for energy in ENERGIES_MEV:
            row = rows.get(float(energy))
            if row is None:
                raise SystemExit("PSTAR %s (%s) has no row at %s MeV" % (number, name, energy))
            picked.append(row)
        materials[key] = {
            "pstar_number": number,
            "name": name,
            "stopping_mev_cm2_g": [row[2] for row in picked],
            "csda_range_g_cm2": [row[3] for row in picked],
        }
    return {
        "about": "NIST PSTAR proton total stopping power and CSDA range (GEN.453). Built by scripts/build_proton_range.py.",
        "source": "NIST Standard Reference Database 124 PSTAR (ICRU Report 49), https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html",
        "retrieved": RETRIEVED,
        "state": "verified",
        "limit": ("Protons to 10 GeV; the nuclear loss (an interaction length of about 90 g/cm2 in air) is not in a CSDA range. "
                  "Above 10 GeV the loader continues the last segment's power law."),
        "energies_mev": list(ENERGIES_MEV),
        "materials": materials,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1], formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--raw", help="a directory of raw PSTAR pages (saved there when missing)")
    parser.add_argument("--check", action="store_true", help="fail if the committed file differs")
    args = parser.parse_args(argv)
    text = json.dumps(build(args.raw), indent=1) + "\n"
    if args.check:
        with open(OUT, encoding="utf-8") as handle:
            if handle.read() != text:
                print("proton_range.json differs from a fresh build", file=sys.stderr)
                return 1
        return 0
    with open(OUT, "w", encoding="utf-8") as handle:
        handle.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
