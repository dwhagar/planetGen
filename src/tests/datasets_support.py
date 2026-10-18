# tests/datasets_support.py

"""Fixtures for the three downloaded datasets of GEN.548, built offline: a small MIST archive, a white-dwarf sequence
archive and an equation-of-state and outer-crust set, each served to the real download step by a fake opener (the sha256 of
the fixture bytes replaces the catalogue's), so no test reaches the network and no downloaded data is committed. Shared by
`test_datasets_readers.py` (the readers), `test_datasets_gate.py` (the refusal) and `test_stellar_composition.py` (GEN.590)."""

import dataclasses
import hashlib
import io
import tarfile
import zipfile

import numpy as np

from planetgen.datasets import catalogue, fetch, isochrones, neutron_star_matter, white_dwarf_models


class _Response(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def _opener(files):
    return lambda request, timeout=None: _Response(files[request.full_url])


def reset_caches(cache_folder, monkeypatch):
    """Points the cache folder at `cache_folder` and forgets everything the readers and the hashing remembered."""
    monkeypatch.setenv(fetch.CACHE_DIR_ENV, str(cache_folder))
    fetch._VERIFIED.clear()
    for loader in (isochrones._load, white_dwarf_models._load, neutron_star_matter._load_compositions,
                   neutron_star_matter._load_eos):
        loader.cache_clear()


def install(monkeypatch, module, dataset, payloads):
    """Points `module.DATASET` at a copy whose pieces are the fixture bytes, and fetches it."""
    pieces = tuple(dataclasses.replace(piece, url="https://example.org/" + piece.name, size=len(payloads[piece.name]),
                                       sha256=hashlib.sha256(payloads[piece.name]).hexdigest())
                   for piece in dataset.pieces)
    fixture = dataclasses.replace(dataset, pieces=pieces)
    monkeypatch.setattr(module, "DATASET", fixture)
    result = fetch.fetch(fixture, opener=_opener({piece.url: payloads[piece.name] for piece in pieces}))
    assert result.status == "fetched", result.reasons
    return fixture


# ---- MIST ----

def _iso_text(fe_h, y=0.27, z=0.014):
    """A MIST-format file: 3 ages; main-sequence masses 0.2 to the turn-off mass, then two evolved rows."""
    lines = ["# MIST version number  = 1.2", "# ---", "#  Yinit        Zinit   [Fe/H]   [a/Fe]  v/vcrit",
             f"# {y}  {z:.5E}     {fe_h:.2f}     0.00     0.40", "# ---", "# number of isochrones =   3"]
    for log_age in (8.0, 9.0, 10.0):
        turn_off = 10 ** ((10 - log_age) / 2.5)
        lines += ["# number of EEPs, cols =   9   25", "#   1 2 3", "# EEP log10_isochrone_age_yr initial_mass"]
        for eep, mass in enumerate(np.linspace(0.2, turn_off, 40)):
            lines.append(_iso_row(eep, log_age, mass, fe_h, 0))
        for extra in range(2):
            lines.append(_iso_row(500 + extra, log_age, turn_off * (1 + 0.001 * extra), fe_h, 2))
        lines.append("")
    return "\n".join(lines).encode()


def _iso_row(eep, log_age, mass, fe_h, phase):
    columns = [0.0] * 25
    columns[0], columns[1], columns[2] = eep, log_age, mass
    columns[7] = 4 * np.log10(mass) + 0.1 * fe_h
    columns[10] = 3.7 + 0.5 * np.log10(mass)
    columns[11] = 0.8 * np.log10(mass)
    columns[24] = phase
    return " ".join(f"{value:.16E}" for value in columns)


def _iso_archive(fe_hs=(-1.0, 0.0)):
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:xz") as archive:
        for fe_h in fe_hs:
            data = _iso_text(fe_h)
            info = tarfile.TarInfo(f"MIST_fixture/MIST_v1.2_feh_{'m' if fe_h < 0 else 'p'}{abs(fe_h):.2f}_basic.iso")
            info.size = len(data)
            archive.addfile(info, io.BytesIO(data))
        note = tarfile.TarInfo("MIST_fixture/README")
        note.size = 3
        archive.addfile(note, io.BytesIO(b"hi\n"))
    return buffer.getvalue()


# ---- white dwarfs ----

def _sequence_text(mass, thin, models=6):
    lines = ["=" * 76, "  #Mod      Teff         Log(g)           R            Age            L",
             "           Log(Tc)       Log(Pc)      Log(rhoc)       Mx/M         Log(qx)",
             "             Lnu        Log(H/*)      Log(He/*)     Log(C/*)      Log(O/*)", "=" * 76]
    for i in range(models):
        teff = 100000.0 / (1 + 4 * i)
        radius = 7e8 * (0.6 / mass) ** (1 / 3) * (1 - 0.05 * i / models)
        lines += [f"{i + 1:6d} {teff:14.4f} 8.0 {radius:.6E} {i * 1e9:.6E} 1.0E+30",
                  "        7.0E+00  2.3E+01  6.5E+00  0.0000000000 0.0E+00",
                  f"        1.0E+30 {-10.0 if thin else -4.0:.6f} -2.0 -0.3054 -0.3054"]
    return "\n".join(lines).encode()


def _sequence_archive(masses=(0.5, 0.6, 0.7)):
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w") as archive:
        for mass in masses:
            for envelope in ("thick", "thin"):
                data = _sequence_text(mass, envelope == "thin")
                info = tarfile.TarInfo(f"seq_{round(mass * 100):03d}_{envelope}.txt")
                info.size = len(data)
                archive.addfile(info, io.BytesIO(data))
    return buffer.getvalue()


# ---- neutron stars ----

def _eos_text(shift):
    rows = ["# rho (fm^-3)  e (MeV/fm^3)  p (MeV/fm^3)"]
    for n in np.geomspace(1e-10, 1.0, 12):
        rows.append(f"{n:.6e} {939.0 * n * (1 + shift):.6e} {n ** (4 / 3) * (10 + shift):.6e}")
    return ("\n".join(rows) + "\n").encode()


def _crust_text(sequence):
    header = "Z1,A1,MASS1[MeV c-2],Z2,A2,MASS2[MeV c-2],nu_max,xe,n1max[fm-3],n2min[fm-3],P12[MeV fm-3],gam12,mu12[MeV],abundance,depth"
    rows = []
    for i, (z, a, density) in enumerate(sequence):
        following = sequence[i + 1] if i + 1 < len(sequence) else None
        z2, a2 = (following[0], following[1]) if following else ("-", "-")
        rows.append(f"{z},{a},1.0e+04,{z2},{a2},1.0e+04,-,1.0,{density:.4e},{density * 1.02:.4e},{density * 0.1:.4e},1.0,930.0,1e-8,0.1")
    return (header + "\n" + "\n".join(rows) + "\n").encode()


def _crust_zip():
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr("compositions/modelA/composition_B0_modelA.csv",
                         _crust_text([(26, 56, 5e-9), (28, 62, 1.6e-7), (28, 64, 8e-7), (36, 86, 1e-4)]))
        archive.writestr("compositions/modelB/composition_B0_modelB.csv", _crust_text([(26, 56, 5e-9), (28, 66, 1e-4)]))
        archive.writestr("compositions/modelA/composition_B1E+14_modelA.csv", b"magnetised, not read\n")
    return buffer.getvalue()


def install_isochrones(monkeypatch):
    return install(monkeypatch, isochrones, isochrones.DATASET, {isochrones.ARCHIVE: _iso_archive()})


def install_white_dwarfs(monkeypatch):
    return install(monkeypatch, white_dwarf_models, white_dwarf_models.DATASET,
                   {white_dwarf_models.ARCHIVE: _sequence_archive()})


def install_neutron_stars(monkeypatch):
    payloads = {neutron_star_matter.COMPOSITIONS: _crust_zip()}
    payloads.update({f"EOS{n}.txt": _eos_text(0.01 * n) for n in range(1, neutron_star_matter.EOS_COUNT + 1)})
    return install(monkeypatch, neutron_star_matter, neutron_star_matter.DATASET, payloads)


def install_all(monkeypatch):
    """All three datasets, on the fixtures, in the cache folder `reset_caches` set; the catalogue lists the fixtures' copies."""
    installed = [install_isochrones(monkeypatch), install_white_dwarfs(monkeypatch), install_neutron_stars(monkeypatch)]
    monkeypatch.setattr(catalogue, "ALL", tuple(installed))
    return installed
