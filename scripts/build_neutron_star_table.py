#!/usr/bin/env python3
# scripts/build_neutron_star_table.py

"""
Builds `src/planetgen/physics/data/neutron_star_curves.json` (GEN.968, wired by GEN.969): the mass-radius curve and the crust of
a neutron star for each of the 21 equations of state of `datasets.neutron_star_matter`, from the TOV integration of
`physics.neutron_star_structure`.

    python scripts/build_neutron_star_table.py            # write the file
    python scripts/build_neutron_star_table.py --check    # integrate again and fail if the committed file differs

Needs the equation-of-state tables in the cache folder (`python3 -m planetgen.cli.datasets fetch neutron_star_matter`) and about
half a minute. The file is small (a few hundred numbers per equation of state); generation reads it and never integrates, so it
needs neither the download nor scipy's integrator at run time.
"""

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

from planetgen.datasets import neutron_star_matter  # noqa: E402
from planetgen.physics import neutron_star_structure as structure  # noqa: E402

OUT = os.path.join(HERE, "..", "src", "planetgen", "physics", "data", "neutron_star_curves.json")
MASSES_PER_CURVE = 25


def build():
    """The table as a dict, ready to dump."""
    models = {}
    for index in range(1, neutron_star_matter.EOS_COUNT + 1):
        masses, radii, crusts, crust_masses = structure.curve_with_crust(index, MASSES_PER_CURVE)
        models[str(index)] = {
            "maximum_mass_solar": round(float(structure.integrated_maximum_mass_solar(index)), 4),
            "mass_solar": [round(float(x), 4) for x in masses],
            "radius_km": [round(float(x), 3) for x in radii],
            "crust_thickness_km": [round(float(x), 3) for x in crusts],
            "crust_mass_solar": [round(float(x), 5) for x in crust_masses],
        }
    return {
        "source": "TOV integration (physics/neutron_star_structure.py) over the 21 equations of state of [R2834] "
                  "(A&A 689, A242, Zenodo [R2834])",
        "crust_core_density_fm3": structure.CRUST_CORE_DENSITY_FM3,
        "models": models,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--check", action="store_true", help="fail if the committed file differs from a new integration")
    arguments = parser.parse_args()
    text = json.dumps(build(), indent=1) + "\n"
    if arguments.check:
        with open(OUT, encoding="utf-8") as handle:
            if handle.read() != text:
                print(f"{os.path.relpath(OUT)} is out of date: run scripts/build_neutron_star_table.py", file=sys.stderr)
                return 1
        return 0
    with open(OUT, "w", encoding="utf-8") as handle:
        handle.write(text)
    print(f"wrote {os.path.relpath(OUT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
