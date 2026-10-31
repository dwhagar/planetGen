# tests/datasets_support.py

"""Fixtures for the downloaded datasets of GEN.548, GEN.443, GEN.715 and GEN.733, built offline: a small MIST archive, a white-dwarf sequence
archive, an equation-of-state and outer-crust set, an ozone cross-section archive, the infrared tables of `ir_datasets_support` and a
hydrogen-helium table archive, each served to the real download step by a fake opener (the sha256 of
the fixture bytes replaces the catalogue's), so no test reaches the network and no downloaded data is committed. Shared by
`test_datasets_readers.py` (the readers), `test_datasets_gate.py` (the refusal) and `test_stellar_composition.py` (GEN.590)."""

import dataclasses
import hashlib
import io
import os
import tarfile
import zipfile

import numpy as np

from planetgen.datasets import (baraffe_isochrones, catalogue, fetch, harris_clusters, hhe_eos, hitran_cia, ir_ktables, isochrones, mamajek_sequence,
                                neutron_star_matter,
                                noaa_proton_events, ozone_uv, white_dwarf_models)
from tests import ir_datasets_support


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
    isochrones.forget()
    for loader in (isochrones._load, white_dwarf_models._load, neutron_star_matter._load_compositions,
                   neutron_star_matter._load_eos, ozone_uv._load, ir_ktables.load, hitran_cia.load, harris_clusters._load,
                   noaa_proton_events._load, baraffe_isochrones._load, mamajek_sequence._load, mamajek_sequence._edges,
                   mamajek_sequence._luminosity_curve, hhe_eos._load):
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

PRE_EDGE_MASSES = {8.0: 0.9, 9.0: 0.35, 10.0: 0.15}
"""dict: The heaviest contracting star (phase -1) of each fixture isochrone age."""


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
        for step, mass in enumerate((PRE_EDGE_MASSES[log_age] * 0.5, PRE_EDGE_MASSES[log_age])):
            lines.append(_iso_row(5 + step, log_age, mass, fe_h, -1))
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


# ---- ozone ----

def ozone_fixture_sigma(nm, temperature_k):
    """The synthetic ozone cross section of the fixture, cm2: a Hartley-like hump at 255 nm and a falling tail, warmer is larger."""
    return (1.1e-17 * np.exp(-((nm - 255.0) / 22.0) ** 2) + 3e-19 * np.exp(-(nm - 300.0) / 18.0) * (nm > 300.0)
            + 5e-21 * np.exp(-((nm - 600.0) / 60.0) ** 2)) * (1.0 + (temperature_k - 218.0) / 770.0)


def ozone_fixture_chappuis_ratio(nm, temperature_k):
    """The fixture's Chappuis band at `temperature_k` over its 298 K self: colder is weaker, most at the blue edge (as Burkholder and Talukdar)."""
    return 1.0 - (298.0 - temperature_k) / 298.0 * 1.2 * np.exp(-(nm - 407.0) / 60.0)


def _ozone_chappuis_text(temperature_k):
    grid = np.arange(407.0, 762.0 + 1e-9, 1.0)
    values = ozone_fixture_sigma(grid, 295.0) * ozone_fixture_chappuis_ratio(grid, temperature_k)
    return ("\n".join(f"{nm:.1f}    {value:.5e}" for nm, value in zip(grid, values)) + "\n\n").encode()


def _ozone_text(low, high, step, temperature_k):
    grid = np.arange(low, high + 1e-9, step)
    return ("\n".join(f"{nm:.2f}    {ozone_fixture_sigma(nm, temperature_k):.5e}" for nm in grid) + "\n\n").encode()


def _ozone_archive():
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        for key, member in ozone_uv.MEMBERS.items():
            if key == "vis":
                archive.writestr(member, _ozone_text(195.0, 830.0, 1.0, 295.0))
            elif key.startswith("bt"):
                archive.writestr(member, _ozone_chappuis_text(float(key[2:])))
            else:
                archive.writestr(member, _ozone_text(195.0, 345.0, 0.05, float(key)))
        archive.writestr("cross_sections/Ozone/readme.txt", b"not read\n")
    return buffer.getvalue()


# ---- hydrogen-helium equation of state ----

HHE_LOG_T = np.arange(2.0, 8.0 + 1e-9, 0.25)
HHE_LOG_P = np.arange(-9.0, 13.0 + 1e-9, 0.25)
"""ndarray: The fixture tables' axes: log10 T in K and log10 P in GPa, a quarter of a dex apart (the real ones are 0.05)."""
HHE_HYDROGEN_MU = 2.3
HHE_HELIUM_MU = 4.0
HHE_STIFFNESS_G_CM3 = 150.0


def hhe_fixture_density_g_cm3(log_t, log_p, mean_molecular_weight):
    """The fixture fluid's density: an ideal gas of that mean molecular weight, held below `HHE_STIFFNESS_G_CM3` by a harmonic
    cap that makes it incompressible at depth, so a body built from it has a finite size."""
    ideal = 1e9 * 10.0 ** log_p * mean_molecular_weight * 1.66054e-27 / (1.380649e-23 * 10.0 ** log_t) / 1000.0
    return 1.0 / (1.0 / ideal + 1.0 / HHE_STIFFNESS_G_CM3)


def hhe_fixture_log_entropy(log_t, log_p):
    """The fixture's log10 entropy (MJ/kg/K): rising with temperature and falling with pressure, so its adiabats run T ~ P^0.4."""
    return -0.5 + 0.2 * (log_t - 4.0) - 0.08 * log_p


def _hhe_text(mean_molecular_weight):
    rows = []
    for log_t in HHE_LOG_T:
        rows.append(f"#iT= log T= {log_t:.3f}")
        for log_p in HHE_LOG_P:
            log_rho = np.log10(hhe_fixture_density_g_cm3(log_t, log_p, mean_molecular_weight))
            values = [log_t, log_p, log_rho, -0.5, hhe_fixture_log_entropy(log_t, log_p), -1.0, 0.5, 5.0, -0.2, 0.4]
            rows.append(" ".join(f"{value:.6E}" for value in values))
    return ("\n".join(rows) + "\n").encode()


def _hhe_archive():
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as archive:
        for member, mu in ((hhe_eos.MIXTURE_MEMBER, 0.725 * HHE_HYDROGEN_MU + 0.275 * HHE_HELIUM_MU), (hhe_eos.HYDROGEN_MEMBER, HHE_HYDROGEN_MU),
                           (hhe_eos.HELIUM_MEMBER, HHE_HELIUM_MU)):
            data = _hhe_text(mu)
            info = tarfile.TarInfo(f"DirEOS2021/{member}")
            info.size = len(data)
            archive.addfile(info, io.BytesIO(data))
        note = tarfile.TarInfo("DirEOS2021/README")
        note.size = 3
        archive.addfile(note, io.BytesIO(b"hi\n"))
    return buffer.getvalue()


def install_hhe_eos(monkeypatch):
    return install(monkeypatch, hhe_eos, hhe_eos.DATASET, {hhe_eos.ARCHIVE: _hhe_archive()})


def install_ozone(monkeypatch):
    return install(monkeypatch, ozone_uv, ozone_uv.DATASET, {ozone_uv.ARCHIVE: _ozone_archive()})


def install_isochrones(monkeypatch):
    return install(monkeypatch, isochrones, isochrones.DATASET, {isochrones.ARCHIVE: _iso_archive()})


def install_white_dwarfs(monkeypatch):
    return install(monkeypatch, white_dwarf_models, white_dwarf_models.DATASET,
                   {white_dwarf_models.ARCHIVE: _sequence_archive()})


def install_neutron_stars(monkeypatch):
    payloads = {neutron_star_matter.COMPOSITIONS: _crust_zip()}
    payloads.update({f"EOS{n}.txt": _eos_text(0.01 * n) for n in range(1, neutron_star_matter.EOS_COUNT + 1)})
    return install(monkeypatch, neutron_star_matter, neutron_star_matter.DATASET, payloads)


# ---- Harris, NOAA proton events, Baraffe ----

HARRIS_TEXT = """\
            Part I:  Identifications and Positional Data

   ID        Name           RA   (2000)   DEC         L       B     R_Sun  R_gc    X     Y     Z

 NGC 104    47 Tuc       00 24 05.67  -72 04 52.6   305.89  -44.89    4.5   7.4   1.9  -2.6  -3.1
 NGC 362                 01 03 14.26  -70 50 55.6   301.53  -46.25    8.6   9.4   3.1  -5.1  -6.2
 Whiting 1               02 02 57     -03 15 10     161.22  -60.76   30.1  34.5 -13.9   4.7 -26.3

_________________________________________________________________________________________________


                    Part II:  Metallicity and Photometry

   ID       [Fe/H] wt  E(B-V) V_HB (m-M)V V_t   M_V,t   U-B   B-V   V-R   V-I  spt   ellip

 NGC 104     -0.72 10   0.04 14.06 13.37  3.95  -9.42   0.37  0.88  0.53  1.14  G4    0.09
 NGC 362     -1.26  5   0.05 15.44 14.83  6.40  -8.43   0.16  0.77  0.49  1.01  F9    0.01
 Whiting 1   -0.70  1   0.03 18.18 17.49 15.03  -2.46                         

_________________________________________________________________________________________________


            Part III:  Velocities and Structural Parameters

    ID         v_r   +/-    v_LSR    sig_v  +/-    c        r_c   r_h    mu_V   rho_0 lg(tc) lg(th)

 NGC 104      -18.0   0.1   -26.7    11.0   0.3   2.07      0.36  3.17   14.38   4.88   7.84  9.55
 NGC 362      223.5   0.5   214.0     6.4   0.3   1.76  c:  0.18  0.82   14.80   4.74   7.76  8.93
 Whiting 1   -130.6   1.8  -138.1                 0.55      0.25  0.22                        7.93

__________________________________________________________________________________________________
"""

NOAA_TEXT = """\
<table border = all>
<tr><th>Start<br>(Day/UT)<th>Maximum<br>(Day/UT)<th>Proton Flux<br>(pfu @ <br>>10 MeV)<th>Year<th>CME<th>Maximum<br>(Day/UT)
<th>Importance<br>(X ray/Opt)<th>Location<th>NOAA SEC<br>Region No.

<tr><td><td><td><td><strong>1976</strong><td><td><td><td><td>

<tr><td>Apr 30/2120<td>May 01/1700<td align = right>    12<td><td><td>Apr 30/2114<td>X2/2B
    <td>S09W47<td align = right> 700

<tr><td><td><td><td><strong>1989</strong><td><td><td><td><td>

<tr><td>Dec 31/2310<td>Jan 01/0400<td align = right>   100<td><td><td>Dec 31/2200<td>M1
    <td>N10W10<td align = right> 1000

<tr><td><td><td><td><strong>2017</strong><td><td><td><td><td>

<tr><td>Aug 10/0115<td>Aug/10 0435<td align = right>    57<td><td><td>Aug 07/1916<td>M4/2B
    <td>S10E2<td align = right>3257

<tr><td>Jun 18/2050<td>June 19/0450<td align = right>    24<td><td>Halo /17 2230<td>Jun 17/2255<td>M6
    <td>S08E61<td align = right>10386

<tr><td>Sep 05/0040<td>Sep&nbsp;08/0035<td align = right>   844<td><td><td>04/2033<td>M5
    <td>S11W16<td align = right>12673
</table>
"""

BARAFFE_TEXT = """\
BHAC15 models

!  t (Gyr) =   0.0010
!-----------------------------------------------------------------------------------------------
! M/Ms  Teff     L/Ls   g    R/Rs   Li/Li0    Mj     Mh     Mk
!-----------------------------------------------------------------------------------------------
 0.010  2400.   -2.40  3.40  0.400    1.000   8.80   8.26   7.85
 0.100  2900.   -1.00  3.40  1.100    1.000   5.60   5.02   4.77
 1.000  4400.    0.60  3.40  3.100    1.000   2.01   1.35   1.22
!-----------------------------------------------------------------------------------------------

!  t (Gyr) =   1.0000
!-----------------------------------------------------------------------------------------------
! M/Ms  Teff     L/Ls   g    R/Rs   Li/Li0    Mj     Mh     Mk
!-----------------------------------------------------------------------------------------------
 0.050  2600.   -3.50  5.00  0.100    0.000   9.00   8.50   8.10
 0.100  2900.   -3.00  5.10  0.130    0.000   8.00   7.50   7.10
 1.000  5700.    0.00  4.40  0.900    0.000   3.00   2.50   2.40
!-----------------------------------------------------------------------------------------------
"""


def install_harris(monkeypatch):
    return install(monkeypatch, harris_clusters, harris_clusters.DATASET, {harris_clusters.FILE: HARRIS_TEXT.encode()})


def install_noaa(monkeypatch):
    return install(monkeypatch, noaa_proton_events, noaa_proton_events.DATASET, {noaa_proton_events.FILE: NOAA_TEXT.encode()})


def install_baraffe(monkeypatch):
    return install(monkeypatch, baraffe_isochrones, baraffe_isochrones.DATASET, {baraffe_isochrones.FILE: BARAFFE_TEXT.encode()})


def _mamajek_text():
    """67 whole-number types O3V to M9V, each 600 K cooler than the one before it down from 44,900 K (not the real table)."""
    types = [f"{letter}{sub}V" for letter, first in zip("OBAFGKM", (3, 0, 0, 0, 0, 0, 0)) for sub in range(first, 10)]
    lines = ["#SpT   Teff  logT   BCv    logL   Mbol R_Rsun   Mv    B-V"]
    for i, name in enumerate(types):
        kelvin = 44900 - 600 * i
        lines.append(f"{name:<6} {kelvin} 4.0 -1.0 {2.0 - 0.05 * i:.2f} 0.0 1.0 1.0 0.0")
    return "\n".join(lines) + "\n"


def install_mamajek(monkeypatch):
    return install(monkeypatch, mamajek_sequence, mamajek_sequence.DATASET, {mamajek_sequence.FILE: _mamajek_text().encode()})


def install_all(monkeypatch):
    """All the datasets, on the fixtures, in the cache folder `reset_caches` set; the catalogue lists the fixtures' copies."""
    installed = [install_isochrones(monkeypatch), install_white_dwarfs(monkeypatch), install_neutron_stars(monkeypatch),
                 install_ozone(monkeypatch), *ir_datasets_support.install(monkeypatch, os.environ[fetch.CACHE_DIR_ENV]),
                 install_harris(monkeypatch), install_noaa(monkeypatch), install_baraffe(monkeypatch), install_mamajek(monkeypatch), install_hhe_eos(monkeypatch)]
    monkeypatch.setattr(catalogue, "ALL", tuple(installed))
    return installed
