#!/usr/bin/env python3
"""
Builds `src/planetgen/physics/data/ozone_chappuis.json`, the ozone cross section of the Chappuis band on the colour model's grid
(380 to 780 nm at 5 nm; GEN.418, `planetgen/physics/atmosphere_colour.py`).

The values are cell means of the 295 K, 1 nm table of the Brion, Daumont, Malicet and Jenouvrier group in the MPI-Mainz UV/VIS Spectral
Atlas snapshot of 2013-06-24 ([R2302, R2303], CC BY 4.0), read through the downloaded dataset `ozone_uv` (run
`python3 -m planetgen.cli.datasets fetch ozone_uv` first), and (GEN.748) the ratio of the band at 220, 240, 260, 280 and 298 K to its
298 K value, cell means of the five tables of [R0599] in the same atlas. Only these derived tables are kept in the
repository, so a page that shows a planet's colour needs no download.

Usage: python3 scripts/build_ozone_chappuis.py
"""

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from planetgen.datasets import fetch, ozone_uv  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src", "planetgen", "physics", "data", "ozone_chappuis.json")
EDGES_NM = np.arange(377.5, 782.6, 5.0)


def main():
    located = fetch.locate(ozone_uv.DATASET)
    if located.missing:
        raise SystemExit("the ozone_uv dataset is not downloaded: run python3 -m planetgen.cli.datasets fetch ozone_uv")
    sigma = ozone_uv.cross_section_cm2(located, EDGES_NM, ozone_uv.CHAPPUIS_TEMPERATURES_K[-1], np.zeros(len(EDGES_NM) - 1))  # the 1 nm table, ratio 1
    centres = 0.5 * (EDGES_NM[:-1] + EDGES_NM[1:])
    ratios = [ozone_uv.chappuis_ratio(located, EDGES_NM, temperature) for temperature in ozone_uv.CHAPPUIS_TEMPERATURES_K]
    document = {
        "source": "the Brion team's 1992 to 1998 measurements, 295 K, 1 nm, MPI-Mainz UV/VIS Spectral Atlas snapshot 2013-06-24, Zenodo record 6951, CC BY 4.0 "
                  "([R2303], Earth Syst. Sci. Data 5, 365); cell means on 5 nm cells",
        "temperature_k": 295.0,
        "wavelength_nm": [float(v) for v in centres],
        "sigma_cm2": [float(f"{v:.6g}") for v in sigma],
        "ratio_source": "[R0599], Geophys. Res. Lett. 21, 581, 407 to 762 nm, MPI-Mainz UV/VIS Spectral Atlas snapshot "
                        "2013-06-24, Zenodo record 6951, CC BY 4.0 [R2302, R2303]; cell means on 5 nm cells over the 298 K table, 1 outside "
                        "407 to 762 nm",
        "ratio_temperatures_k": list(ozone_uv.CHAPPUIS_TEMPERATURES_K),
        "ratio": [[float(f"{v:.5g}") for v in row] for row in ratios],
    }
    with open(OUT, "w", encoding="utf-8") as handle:
        json.dump(document, handle, indent=1)
        handle.write("\n")
    print(f"wrote {os.path.relpath(OUT)}: peak {max(document['sigma_cm2']):.3e} cm2 at {centres[int(np.argmax(sigma))]:.1f} nm")


if __name__ == "__main__":
    main()
