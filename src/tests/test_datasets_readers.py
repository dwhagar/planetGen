# tests/test_datasets_readers.py

"""The three datasets of GEN.548 read through their functions, on fixtures the tests build offline: a small MIST
archive, a white-dwarf sequence archive and an equation-of-state and outer-crust set. No test reaches the network."""

import dataclasses
import hashlib
import io
import tarfile
import zipfile

import numpy as np
import pytest

from planetgen.datasets import catalogue, fetch, isochrones, neutron_star_matter, white_dwarf_models
from planetgen.physics import constants


class _Response(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def _opener(files):
    return lambda request, timeout=None: _Response(files[request.full_url])


@pytest.fixture(autouse=True)
def cache(tmp_path, monkeypatch):
    monkeypatch.setenv(fetch.CACHE_DIR_ENV, str(tmp_path / "cache"))
    fetch._VERIFIED.clear()
    for loader in (isochrones._load, white_dwarf_models._load, neutron_star_matter._load_compositions,
                   neutron_star_matter._load_eos):
        loader.cache_clear()
    return tmp_path / "cache"


def _install(monkeypatch, module, dataset, payloads):
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


@pytest.fixture
def mist(monkeypatch):
    return _install(monkeypatch, isochrones, isochrones.DATASET, {isochrones.ARCHIVE: _iso_archive()})


def test_the_derived_file_keeps_only_main_sequence_rows(mist, cache):
    assert not (cache / "isochrones" / isochrones.ARCHIVE).exists()
    with np.load(cache / "isochrones" / isochrones.DERIVED) as data:
        assert list(data["fe_h"]) == [-1.0, 0.0] and data["rows_0"].shape == (120, 5)
        assert data["rows_0"][:, 1].max() <= 10 ** ((10 - 8.0) / 2.5) + 1e-3


def test_a_main_sequence_point_interpolates_in_mass_age_and_metallicity(mist):
    log_l, log_teff, log_r = isochrones.main_sequence_point(0.5, 1.0, 0.0)
    assert log_l == pytest.approx(4 * np.log10(0.5), abs=0.02)
    assert log_teff == pytest.approx(3.7 + 0.5 * np.log10(0.5), abs=0.005)
    assert log_r == pytest.approx(0.8 * np.log10(0.5), abs=0.01)
    assert isochrones.main_sequence_point(0.5, 1.0, -0.5)[0] == pytest.approx(4 * np.log10(0.5) - 0.05, abs=0.02)
    between = isochrones.main_sequence_point(0.5, 0.5, 0.0)
    assert between is not None and between[0] == pytest.approx(4 * np.log10(0.5), abs=0.02)


def test_a_star_that_left_the_main_sequence_or_is_outside_the_ages_gets_none(mist):
    assert isochrones.main_sequence_point(2.0, 5.0, 0.0) is None     # turn-off at 5 Gy is 1.2 Msun
    assert isochrones.main_sequence_point(0.5, 1e-5, 0.0) is None    # younger than the grid
    assert isochrones.main_sequence_point(0.5, 200.0, 0.0) is None   # older than the grid
    assert isochrones.main_sequence_point(0.1, 1.0, 0.0) is None     # lighter than the grid's lowest mass


def test_the_lifetime_follows_the_turn_off_and_metallicity_is_clamped(mist):
    assert isochrones.main_sequence_lifetime_gy(10 ** (1 / 2.5), 0.0) == pytest.approx(1.0, rel=0.02)  # a grid age
    between = isochrones.main_sequence_lifetime_gy(2.0, 0.0)
    assert 1.0 < between < 10.0
    assert isochrones.main_sequence_lifetime_gy(0.5, 0.0) is None
    assert isochrones.main_sequence_point(0.5, 1.0, 3.0) == isochrones.main_sequence_point(0.5, 1.0, 0.0)
    assert isochrones.metallicities() == (-1.0, 0.0)
    assert isochrones.abundances(-0.5) == pytest.approx((0.27, 0.014))


def test_the_parser_refuses_a_file_that_is_not_a_mist_isochrone():
    with pytest.raises(ValueError):
        isochrones.parse_isochrone_file([b"# nothing\n", b"1 2 3\n"])


def test_the_isochrone_readers_say_so_when_the_archive_was_never_downloaded(cache):
    with pytest.raises(fetch.StandIn, match=r"stand-in \(SC1\)"):
        isochrones.main_sequence_point(1.0, 1.0, 0.0)


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


@pytest.fixture
def dwarfs(monkeypatch):
    return _install(monkeypatch, white_dwarf_models, white_dwarf_models.DATASET,
                    {white_dwarf_models.ARCHIVE: _sequence_archive()})


def test_the_sequence_parser_reads_sixteen_numbers_per_model(dwarfs):
    first = white_dwarf_models.sequence(0.6)
    assert len(first.teff_k) == 6 and first.teff_k[0] == pytest.approx(100000.0)
    assert first.radius_cm[0] == pytest.approx(7e8)
    assert first.log_mass_fractions["hydrogen"][0] == pytest.approx(-4.0)
    assert first.log_mass_fractions["oxygen"][0] == pytest.approx(-0.3054)


def test_the_radius_interpolates_in_temperature_and_mass(dwarfs):
    one_solar_radius_cm = constants.SOLAR_RADIUS_M * 100.0
    exact = white_dwarf_models.radius_solar(0.6, 100000.0)
    assert exact == pytest.approx(7e8 / one_solar_radius_cm)
    middle = white_dwarf_models.radius_solar(0.55, 100000.0)
    assert white_dwarf_models.radius_solar(0.6, 100000.0) < middle < white_dwarf_models.radius_solar(0.5, 100000.0)
    assert white_dwarf_models.radius_solar(0.6, 1e9) == pytest.approx(exact)  # hotter than the sequence: its hot end
    assert white_dwarf_models.radius_solar(0.6, 1.0) < exact


def test_masses_outside_the_grid_get_none_and_the_masses_are_listed(dwarfs):
    assert white_dwarf_models.masses_sol() == (0.5, 0.6, 0.7)
    assert white_dwarf_models.radius_solar(0.4, 1e4) is None and white_dwarf_models.radius_solar(0.8, 1e4) is None


def test_the_layers_are_the_cases_and_the_core_is_split_evenly(dwarfs):
    thick = white_dwarf_models.layer_mass_fractions(0.6)
    thin = white_dwarf_models.layer_mass_fractions(0.6, "thin")
    assert thick["hydrogen"] == pytest.approx(1e-4) and thin["hydrogen"] == pytest.approx(1e-10)
    assert thick["helium"] == pytest.approx(1e-2) and thick["carbon"] == thick["oxygen"]
    assert sum(thick.values()) == pytest.approx(1.0, abs=0.02)


def test_a_sequence_file_with_a_broken_model_is_refused():
    with pytest.raises(ValueError):
        white_dwarf_models.parse_sequence("\n".join(["x"] * 5 + ["1.0 2.0 3.0"]), 0.6, "thick")


def test_the_white_dwarf_readers_say_so_when_nothing_was_downloaded(cache):
    with pytest.raises(fetch.StandIn, match=r"stand-in \(SC2\)"):
        white_dwarf_models.radius_solar(0.6, 1e4)


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


@pytest.fixture
def stars(monkeypatch):
    payloads = {neutron_star_matter.COMPOSITIONS: _crust_zip()}
    payloads.update({f"EOS{n}.txt": _eos_text(0.01 * n) for n in range(1, neutron_star_matter.EOS_COUNT + 1)})
    return _install(monkeypatch, neutron_star_matter, neutron_star_matter.DATASET, payloads)


def test_all_equations_of_state_are_read_and_ascend_in_density(stars):
    tables = neutron_star_matter.equations_of_state()
    assert [table.index for table in tables] == list(range(1, 22))
    assert all(np.all(np.diff(table.baryon_density_fm3) > 0) for table in tables)
    assert tables[0].pressure_mev_fm3[-1] < tables[-1].pressure_mev_fm3[-1]


def test_an_equation_of_state_with_falling_density_is_refused():
    with pytest.raises(ValueError):
        neutron_star_matter.parse_equation_of_state(1, "# h\n2 1 1\n1 1 1\n")


def test_only_the_unmagnetised_compositions_are_read(stars):
    assert neutron_star_matter.crust_models() == ("modelA", "modelB")
    layers = neutron_star_matter.outer_crust("modelA")
    assert [(layer.z, layer.a) for layer in layers] == [(26, 56), (28, 62), (28, 64), (36, 86)]


def test_the_crust_nucleus_follows_the_density_and_ends_at_drip(stars):
    grams_per_cm3 = lambda n_fm3: n_fm3 / 1e-39 * 1.66053906660e-24
    assert neutron_star_matter.crust_nucleus("modelA", grams_per_cm3(1e-9)) == (26, 56)
    assert neutron_star_matter.crust_nucleus("modelA", grams_per_cm3(1e-7)) == (28, 62)
    assert neutron_star_matter.crust_nucleus("modelA", grams_per_cm3(5e-7)) == (28, 64)
    assert neutron_star_matter.crust_nucleus("modelB", grams_per_cm3(5e-5)) == (28, 66)
    assert neutron_star_matter.crust_nucleus("modelA", grams_per_cm3(1e-3)) is None
    assert neutron_star_matter.baryon_density_fm3(grams_per_cm3(1e-9)) == pytest.approx(1e-9)


def test_the_neutron_star_readers_say_so_when_nothing_was_downloaded(cache):
    with pytest.raises(fetch.StandIn, match=r"stand-in \(SC3\)"):
        neutron_star_matter.crust_models()


# ---- the catalogue ----

def test_the_catalogue_lists_the_three_datasets_with_distinct_fills():
    assert [dataset.key for dataset in catalogue.ALL] == ["isochrones", "white_dwarf_models", "neutron_star_matter"]
    assert [dataset.fill for dataset in catalogue.ALL] == ["SC1", "SC2", "SC3"]
    assert len(neutron_star_matter.DATASET.pieces) == 1 + neutron_star_matter.EOS_COUNT
    assert all(piece.url.startswith("https://") for dataset in catalogue.ALL for piece in dataset.pieces)
