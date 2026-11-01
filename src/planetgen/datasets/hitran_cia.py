# planetgen/datasets/hitran_cia.py

"""
Collision-induced absorption of the gas pairs the greenhouse model needs (GEN.715, AD11): the binary absorption coefficient
`k(nu, T)` in cm5 per molecule2, reduced to the band means of the infrared opacity.

Source: the HITRAN collision-induced absorption files ([R2267]; the
2011, 2021 and 2024 releases named below), from hitran.org/data/CIA/main:

- N2-N2 (2021), O2-O2 (2024), N2-O2 (2024), CO2-CO2 (2024), N2-CH4 (2024), N2-H2 (2024), CO2-CH4 (2024), CO2-H2 (2024),
- H2-H2 (2011), H2-CH4 (2011, normal hydrogen), H2-He (2011).

The files are transient pieces: setup downloads each once, checks its sha256, writes the mean of `k` over each of the 130 bands
of `EDGES_CM1` for every temperature of the file into `ir_cia.npz`, and removes them. A file is a set of blocks, each
one temperature over one wavenumber range; the ranges differ with the temperature (the carbon dioxide file has its far-infrared block at 200,
224, 250 K and so on, its 1150 to 1850 cm-1 block at others). For a band the blocks that cover at least half of it are used, a temperature's
value is the mean of `k` over the part the block covers, a band no block covers is zero, and the values between and beyond the temperatures
are held at the nearest measured one (the reader interpolates linearly in temperature between them).

Limits (stated; registry row IR6 in docs/design/research/gap-filling-plan.md): a band measured at a few temperatures only
(oxygen in the near infrared) is flat in temperature; pairs the files do not hold (CH4-CH4, the noble gases' pairs other than H2-He, every
pair with a water molecule) are zero, the water continuum being the separate `line_opacity.WATER_CONTINUUM`.
"""

import functools
import io

import numpy as np

from planetgen.datasets import fetch, ir_ktables

DERIVED = "ir_cia.npz"
SOURCE = "HITRAN collision-induced absorption [R2267]"
URL = "https://hitran.org/data/CIA/main/"

FILES = {
    "N2-N2": ("N2-N2_2021.cia", "f0525f75c667802e8dbbfd7054a0a7e07863d8a4ffdb449fe80666bc658c849d", 5692968),
    "O2-O2": ("O2-O2_2024.cia", "0146acfe5c1414b8deaa56be9a63159154479d5804d1615855888a2e1efd4961", 4203198),
    "N2-O2": ("O2-N2_2024.cia", "90851fbfae226d7ecda93efc0de16105f72966a532adc7f1bc234cb13833799f", 1775796),
    "CO2-CO2": ("CO2-CO2_2024.cia", "9db0856cb008cb147728e3a945f818c60a02370442d492c38d7f8884c9cfadca", 2115361),
    "CH4-N2": ("N2-CH4_2024.cia", "7056ee0bc4a38a4a4d2b215b40312ba5be6029dd6098bf1b1c76d5522eb12523", 6532284),
    "H2-N2": ("N2-H2_2024.cia", "986d73551da83810d4729cfaa4f49a7256eea83c6fc95ec51b4ffc1bb80badb1", 449814),
    "CH4-CO2": ("CO2-CH4_2024.cia", "e6fde5acdbe7586f08c98adfc5ad0c694051cfe20ffad654ea62bbce00a69e17", 2442276),
    "CO2-H2": ("CO2-H2_2024.cia", "f547b53f89108e064a8bb03316fad4588ead058a892673035c127f2e8a7f50d7", 43812),
    "H2-H2": ("H2-H2_2011.cia", "14f56ff889a474402136be721fcc286921a1b47ce5f7c3dff90af2d4c74b20e5", 24824179),
    "CH4-H2": ("H2-CH4_norm_2011.cia", "7c32a42d1491af41de9891c32a9c729c36e90814da090c045382610b0c13d833", 435290),
    "H2-He": ("H2-He_2011.cia", "4f0eb9cd69a1c383f53a1431495bae0c01a30a41cc1b8433c2d726763ff45431", 146854122),
}
"""dict: A pair (its two formulas sorted by name, joined by `-`) to (the file at hitran.org, sha256, size in bytes)."""

EDGES_CM1 = np.concatenate([[0.0], ir_ktables.BAND_EDGES_CM1])
"""ndarray: The 131 edges, cm-1, of the 130 bands: the k-tables' 129 and one from 0 to 200 cm-1 under them."""

COVERAGE = 0.5
"""float: The share of a band a block must cover to count for the band."""


def pair_key(first, second):
    """The key of a pair of formulas in `FILES`, `ir_cia.npz` and the readers: the two formulas sorted, joined by `-`."""
    return "-".join(sorted((first, second)))


def parse(text):
    """
    The blocks of a HITRAN CIA file: a list of `(temperature_k, array of (wavenumber cm-1, k cm5/molecule2))`.

    Raises:
        ValueError: A header or a data line is not the format.
    """
    lines = text.splitlines()
    blocks = []
    index = 0
    while index < len(lines):
        header = lines[index]
        if not header.strip():
            index += 1
            continue
        try:
            fields = header[20:].split()
            count, temperature = int(fields[2]), float(fields[3])
            data = np.array([[float(value) for value in line.split()[:2]] for line in lines[index + 1:index + 1 + count]])
        except (IndexError, ValueError) as exc:
            raise ValueError(f"line {index + 1}: not a CIA block header or its data ({exc})") from exc
        if data.shape != (count, 2) or np.any(np.diff(data[:, 0]) <= 0):
            raise ValueError(f"line {index + 1}: a CIA block must hold {count} rising wavenumbers")
        data[:, 1] = np.maximum(data[:, 1], 0.0)        # the files hold noise a little under zero in the wings
        blocks.append((temperature, data))
        index += 1 + count
    if not blocks:
        raise ValueError("the CIA file holds no block")
    return blocks


def band_means(blocks, edges_cm1):
    """
    The mean of `k` over each band of `edges_cm1` at every temperature of the blocks.

    Returns:
        tuple: `(temperature_k, k)` with `k` of shape (temperature, band) in cm5 per molecule2, the temperatures the union of the
        blocks' (sorted, unique), each band filled between and beyond the temperatures a block covers it at.
    """
    edges = np.asarray(edges_cm1, dtype=float)
    width = np.diff(edges)
    temperatures = np.array(sorted({temperature for temperature, _data in blocks}))
    sums = np.zeros((len(temperatures), len(width)))
    counts = np.zeros((len(temperatures), len(width)))
    for temperature, data in blocks:
        row = int(np.searchsorted(temperatures, temperature))
        wavenumber, k = data[:, 0], data[:, 1]
        cumulative = np.concatenate([[0.0], np.cumsum(0.5 * (k[1:] + k[:-1]) * np.diff(wavenumber))])
        low = np.clip(edges[:-1], wavenumber[0], wavenumber[-1])
        high = np.clip(edges[1:], wavenumber[0], wavenumber[-1])
        covered = high - low
        usable = covered >= COVERAGE * width
        area = np.interp(high, wavenumber, cumulative) - np.interp(low, wavenumber, cumulative)
        mean = np.where(covered > 0.0, area / np.where(covered > 0.0, covered, 1.0), 0.0)
        sums[row] += np.where(usable, mean, 0.0)
        counts[row] += usable
    result = np.zeros_like(sums)
    for band in range(len(width)):
        have = counts[:, band] > 0
        if have.any():
            result[:, band] = np.interp(temperatures, temperatures[have], sums[have, band] / counts[have, band])
    return temperatures, result


def derive(folder):
    """Writes `ir_cia.npz` from the verified files in `folder`: `temperature_<pair>` and `k_<pair>` for each pair of `FILES`."""
    arrays = {}
    for pair, (name, _sha, _size) in FILES.items():
        with open(f"{folder}/{name}", encoding="ascii", errors="replace") as handle:
            blocks = parse(handle.read())
        arrays[f"temperature_{pair}"], arrays[f"k_{pair}"] = band_means(blocks, EDGES_CM1)
    buffer = io.BytesIO()
    np.savez_compressed(buffer, **arrays)
    return {DERIVED: buffer.getvalue()}


DATASET = fetch.Dataset(
    key="ir_cia",
    title="Collision-induced absorption of N2, O2, CO2, H2, He and CH4 pairs (HITRAN CIA, 2011 to 2024 releases)",
    source="[R2267], HITRAN collision-induced absorption files from hitran.org",
    pieces=tuple(fetch.Piece(name, URL + name, sha, size, transient=True) for name, sha, size in FILES.values()),
    replaces="the coarse pair coefficients of the grey infrared opacity, read at six temperatures from a part-03 table",
    fill="AD11",
    fallback="none at generation: the fitted greenhouse table (physics/greenhouse_table.py) carries the answer; the pairs are read to fit and check it",
    derived=(DERIVED,),
    derive=lambda folder: derive(folder),
    derive_version=1,
)


@functools.lru_cache(maxsize=2)
def load(path):
    with np.load(path) as stored:
        return {name: np.array(stored[name]) for name in stored.files}


def tables(located):
    """The band means as a dict of arrays (`temperature_<pair>`, `k_<pair>`) from a `fetch.Located`."""
    return load(located.path(DERIVED))
