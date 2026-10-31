#!/usr/bin/env python3
"""
Builds `src/planetgen/physics/data/acid_cloud_optics.json` (GEN.1020): the infrared optics of sulfuric acid cloud droplets on the 130 bands of
`physics/line_opacity.py`, from the refractive indices of the HITRAN aerosol package (Myhre and others 2003 for 76 and 81 percent acid at 298 K,
400 to 7500 cm-1; Palmer and Williams 1975 for the visible index).

    python scripts/build_acid_cloud_optics.py --indices /path/to/hitran_ri      # the folder of hitran_ri.tar (https://hitran.org/data/Aerosols/)

For each acid strength and each effective radius of `RADII_UM` the droplets are a gamma size distribution of effective variance `VARIANCE`
(0.2, the order of the mode variances of Arney and Kane 2018 Table 2); the Mie efficiencies (`physics/mie.py`) are averaged over
four wavenumbers of each band. Stored per band: the extinction efficiency relative to the 550 nm one of the same droplets (the cloud depth in
the generator is stated at 550 nm), the single-scattering albedo and the asymmetry factor. The raw index files are never committed.
"""

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

from planetgen.physics import line_opacity, mie  # noqa: E402

RADII_UM = (0.3, 0.5, 0.7, 1.0, 1.4, 2.0, 3.0, 3.85)
VARIANCE = 0.2
SAMPLES_PER_BAND = 4
STRENGTHS = {76: ("rnw76t298", 3), 81: ("rnw81t298", 0)}
"""dict: Weight percent to the label and the column of the 298 K set in `myhreh2so4set1.dat`."""
VISIBLE_CM1 = 18000.0
OUTPUT = os.path.join(HERE, "..", "src", "planetgen", "physics", "data", "acid_cloud_optics.json")


def _read_myhre(folder):
    """`(wavenumber, {percent: complex index})` from set 1, ascending wavenumber."""
    with open(os.path.join(folder, "ascii", "myhre_h2so4", "myhreh2so4set1.dat"), encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    real = [line.split() for line in lines[15:15 + 4733]]
    imaginary = [line.split() for line in lines[4750:4750 + 4733]]
    nu = np.array([float(r[0]) for r in real])
    result = {}
    for percent, (_label, column) in STRENGTHS.items():
        n = np.array([float(r[2 + column]) for r in real])
        k = np.array([float(r[2 + column]) for r in imaginary])
        result[percent] = n + 1j * k
    order = np.argsort(nu)
    return nu[order], {percent: index[order] for percent, index in result.items()}


def _visible_index(folder):
    """The real index at 550 nm of 75 and 84.5 percent acid (Palmer and Williams: the imaginary part is zero there), averaged."""
    with open(os.path.join(folder, "ascii", "single_files", "palmer_williams_h2so4.dat"), encoding="utf-8") as handle:
        for line in handle:
            cells = line.split()
            if cells and cells[0] == "18000.":
                return 0.5 * (float(cells[5]) + float(cells[6]))
    raise ValueError("no 18000 cm-1 row in the Palmer and Williams file")


def build(folder):
    nu_data, indices = _read_myhre(folder)
    visible_n = _visible_index(folder)
    edges = line_opacity.EDGES_CM1
    sample_nu = np.concatenate([edges[i] + (edges[i + 1] - edges[i]) * (np.arange(SAMPLES_PER_BAND) + 0.5) / SAMPLES_PER_BAND
                                for i in range(len(edges) - 1)])
    sample_nu = np.maximum(sample_nu, 1.0)
    document = {
        "version": 1,
        "title": "Infrared optics of sulfuric acid cloud droplets on the infrared model's bands",
        "description": "Derived by scripts/build_acid_cloud_optics.py from the HITRAN aerosol refractive indices of Myhre and others 2003 (76 and 81 percent H2SO4 at 298 K, 400 to 7500 cm-1, held at the ends) and Palmer and Williams 1975 (the 550 nm index), with Mie spheres in a gamma size distribution of effective variance 0.2. Per band: extinction efficiency over the 550 nm one of the same droplets, single-scattering albedo, asymmetry factor.",
        "variance": VARIANCE, "radii_um": list(RADII_UM), "weight_percent": list(STRENGTHS), "bands": len(edges) - 1,
        "visible_real_index": visible_n, "optics": {},
    }
    for percent in STRENGTHS:
        index = np.interp(sample_nu, nu_data, indices[percent].real) + 1j * np.interp(sample_nu, nu_data, indices[percent].imag)
        rows = []
        for radius in RADII_UM:
            wavelength = 1.0e4 / sample_nu
            q_ext, q_sca, g = mie.gamma_average(wavelength, index, radius, VARIANCE)
            q_visible = mie.gamma_average(np.array([1.0e4 / VISIBLE_CM1]), np.array([visible_n + 0j]), radius, VARIANCE)[0][0]
            q_ext, q_sca, g = (a.reshape(-1, SAMPLES_PER_BAND) for a in (q_ext, q_sca, g))
            ext = q_ext.mean(axis=1)
            sca = q_sca.mean(axis=1)
            asym = np.where(sca > 0.0, (q_sca * g).mean(axis=1) / np.maximum(sca, 1e-300), 0.0)
            rows.append({"q_visible": round(float(q_visible), 4), "ext_over_visible": [round(float(v), 4) for v in ext / q_visible],
                         "omega": [round(float(v), 4) for v in np.where(ext > 0.0, sca / ext, 0.0)],
                         "g": [round(float(v), 4) for v in asym]})
        document["optics"][str(percent)] = rows
    return document


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("--indices", required=True, help="the folder of hitran_ri.tar (it holds ascii/)")
    args = parser.parse_args(argv)
    document = build(args.indices)
    with open(OUTPUT, "w", encoding="utf-8") as handle:
        json.dump(document, handle, indent=1)
        handle.write("\n")
    print(f"wrote {os.path.relpath(OUTPUT)}")


if __name__ == "__main__":
    main()
