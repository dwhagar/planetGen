# planetgen/physics/neutron_star_crust.py

"""
The matter of the crust of a neutron star by density (GEN.1358)
===============================================================

The crust of a neutron star is a lattice of nuclei in an electron gas whose nuclei change with the density, down to the layer where the matter turns
uniform. [R0722: Tables 2 and 3] (READ) gives the sequence, and `data/neutron_star_crust.json` holds it:

* **The outer crust**, Table 2: a single kind of nucleus at each density, the ground state of cold catalysed matter, from iron-56 at the surface
  through nickel and the nuclei with 50 neutrons (krypton-86, selenium-84, germanium-82, zinc-80, zinc-82) to the nuclei with 82 neutrons (palladium-128
  down to krypton-118). The masses are experimental up to zinc-82 and the BSk8 mass model beyond. Each nucleus holds from the density where it replaces
  the one above up to `max_density_g_cm3`.
* **Neutron drip** at 4.27e11 g/cm3, where the nuclei hold all the neutrons they can and the rest run free.
* **The inner crust**, Table 3 (the Negele and Vautherin calculation): at 11 densities from 4.67e11 to 1.3e14 g/cm3 a cell holds one cluster in a gas of
  free neutrons: zirconium (Z = 40) to 1.5e12, tin (Z = 50) from 6e12 to 3.4e13, zirconium and germanium toward the bottom. The labels count the
  nucleons of the whole cell, most of them the free gas, not the cluster.

* **The nuclear pasta** (GEN.1427), the bottom layer of the crust: where the nuclei crowd (from 0.2 of the normal nuclear density, `PASTA`) the energy of
  the system falls by merging them into rods, slabs and the inverse shapes, so the cells of the Negele and Vautherin table above give way to gnocchi,
  spaghetti, lasagna, then the inverses (cylindrical bubbles, then spherical bubbles of the neutron gas) down to the crust-core boundary, `neutron_star_structure.CRUST_CORE_DENSITY_FM3`
  (0.08 fm^-3, row NS1). The layer is about 100 m thick and may hold half the mass of the crust [R5167] (READ). It is a prediction of calculations, not an
  observation (evidence status: predicted only), so the page says so; there is one layer row, the order of the shapes is the review's, and the densities where one
  shape gives way to the next are not given by it.

The core of neutrons, protons and electrons beyond is not in the table. There is no stable bulk neutron solid outside a star, so the
matter has no row in the materials table.

The sequence is the same for every neutron star (no star-dependent draw), so the module holds data and readers only: the neutron star page shows it
(`describe.build`).
"""

import json
import os
from dataclasses import dataclass

from planetgen.physics import elements, neutron_star_structure
from planetgen.physics.elements import DataError

G_CM3_PER_FM3 = 1.66053906660e-24 * 1.0e39
"""float: One baryon per fm^3 in g/cm3 (the atomic mass unit in g times the 1e39 fm^3 in a cm3), to give the pasta's densities in the table's unit."""

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "neutron_star_crust.json")
"""str: The data file."""


@dataclass(frozen=True)
class OuterNucleus:
    """The nucleus of one layer of the outer crust: `max_density_g_cm3` is the density where the next one replaces it."""

    symbol: str
    mass_number: int
    protons: int
    neutrons: int
    max_density_g_cm3: float
    mass_data: str

    @property
    def name(self):
        """`"iron-56"`."""
        return f"{elements.element(self.symbol).name.lower()}-{self.mass_number}"


@dataclass(frozen=True)
class PastaShape:
    """One shape of the nuclear pasta, `name` its nickname and `shape` what the matter looks like."""

    name: str
    shape: str


@dataclass(frozen=True)
class Pasta:
    """The nuclear pasta layer (GEN.1427): `shapes` from the top of the layer down (rising density), the density range `lower_g_cm3` to `upper_g_cm3`, g/cm3, the
    thickness `thickness_m`, m, the most of the crust's mass it may hold, `crust_mass_share_at_most`, and the evidence status."""

    shapes: tuple
    lower_g_cm3: float
    upper_g_cm3: float
    thickness_m: float
    crust_mass_share_at_most: float
    evidence: str


@dataclass(frozen=True)
class InnerCell:
    """A cell of the inner crust at `density_g_cm3`: a cluster of `protons` protons (the element `symbol`) in a gas of free neutrons; the nucleons
    of the whole cell are `nucleons_in_cell`, `neutrons_in_cell` of them neutrons."""

    density_g_cm3: float
    symbol: str
    nucleons_in_cell: int
    protons: int
    neutrons_in_cell: int

    @property
    def name(self):
        """`"zirconium"`: the element of the cluster."""
        return elements.element(self.symbol).name.lower()


def load(path=DATA_FILE):
    """Reads and checks the file at `path`: `(outer crust, inner crust, neutron drip density, pasta)`; a row that does not add up is refused."""
    try:
        with open(path, encoding="utf-8") as handle:
            document = json.load(handle)
    except (OSError, ValueError) as error:
        raise DataError(f"cannot read the neutron star crust table {path}: {error}") from error

    def number(value, what):
        if not isinstance(value, (int, float)) or isinstance(value, bool) or not value > 0 or value == float("inf"):
            raise DataError(f"{what} must be a positive number, not {value!r}")
        return value

    outer = []
    for raw in document.get("outer_crust") or []:
        symbol = raw.get("symbol")
        protons, neutrons, mass = int(number(raw.get("protons"), "the protons")), int(number(raw.get("neutrons"), "the neutrons")), \
            int(number(raw.get("mass_number"), "the mass number"))
        if protons + neutrons != mass or elements.element(protons).symbol != symbol:
            raise DataError(f"outer crust nucleus {mass}{symbol} does not add up: {protons} protons and {neutrons} neutrons")
        outer.append(OuterNucleus(symbol, mass, protons, neutrons, float(number(raw.get("max_density_g_cm3"), "the density")), raw.get("mass_data", "")))
    inner = []
    for raw in document.get("inner_crust") or []:
        symbol = raw.get("symbol")
        protons, neutrons, nucleons = int(number(raw.get("protons"), "the protons")), int(number(raw.get("neutrons_in_cell"), "the neutrons")), \
            int(number(raw.get("nucleons_in_cell"), "the nucleons"))
        if protons + neutrons != nucleons or elements.element(protons).symbol != symbol:
            raise DataError(f"inner crust cell {nucleons}{symbol} does not add up: {protons} protons and {neutrons} neutrons")
        inner.append(InnerCell(float(number(raw.get("density_g_cm3"), "the density")), symbol, nucleons, protons, neutrons))
    if len(outer) < 10 or len(inner) < 5:
        raise DataError("the neutron star crust table is missing layers")
    for rows, density in ((outer, lambda row: row.max_density_g_cm3), (inner, lambda row: row.density_g_cm3)):
        values = [density(row) for row in rows]
        if values != sorted(values) or len(set(values)) != len(values):
            raise DataError("the densities of a layer table must rise, one row each")
    drip = float(number(document.get("neutron_drip_g_cm3"), "the neutron drip density"))
    if drip != outer[-1].max_density_g_cm3 or inner[0].density_g_cm3 <= drip:
        raise DataError("the outer crust ends at neutron drip and the inner crust starts below it")
    raw = document.get("pasta") or {}
    shapes = tuple(PastaShape(str(item.get("name") or ""), str(item.get("shape") or "")) for item in raw.get("shapes") or [])
    if len(shapes) < 3 or not all(item.name and item.shape for item in shapes):
        raise DataError("the pasta needs its shapes, each with a name and a description")
    saturation = float(number(raw.get("saturation_density_fm3"), "the normal nuclear density"))
    fraction = float(number(raw.get("lower_fraction_of_saturation"), "the pasta's lower density fraction"))
    lower = fraction * saturation * G_CM3_PER_FM3
    if not lower < neutron_star_structure.CRUST_CORE_DENSITY_FM3 * G_CM3_PER_FM3:
        raise DataError("the pasta starts above the crust-core boundary")
    share = float(number(raw.get("crust_mass_share_at_most"), "the pasta's share of the crust mass"))
    if share > 1.0:
        raise DataError("the pasta's share of the crust mass is a fraction")
    pasta = Pasta(shapes, lower, neutron_star_structure.CRUST_CORE_DENSITY_FM3 * G_CM3_PER_FM3, float(number(raw.get("thickness_m"), "the pasta's thickness")), share,
                  str(raw.get("evidence") or ""))
    return tuple(outer), tuple(inner), drip, pasta


OUTER_CRUST, INNER_CRUST, NEUTRON_DRIP_G_CM3, PASTA = load()
"""The nuclei of the outer crust from the surface down (`OuterNucleus`), the cells of the inner crust from the top (`InnerCell`), the density
of neutron drip, g/cm3, and the nuclear pasta layer at the bottom (`Pasta`, GEN.1427)."""

CRUST_BOTTOM_G_CM3 = INNER_CRUST[-1].density_g_cm3
"""float: The density of the deepest tabulated cell, g/cm3 (the table stops short of the crust-core transition)."""


def inner_runs():
    """The inner crust as runs of one element: `[(name, first density, last density)]`, g/cm3, from the top down (a run of one cell has the same
    density twice)."""
    runs = []
    for cell in INNER_CRUST:
        if runs and runs[-1][0] == cell.name:
            runs[-1][2] = cell.density_g_cm3
        else:
            runs.append([cell.name, cell.density_g_cm3, cell.density_g_cm3])
    return [tuple(run) for run in runs]


def sequence():
    """
    The crust from the surface down as `(region, name, upper density g/cm3 or None)` rows: each outer crust nucleus with the density up to which it holds,
    then `("drip", "free neutrons", drip density)`, then each inner crust cell as its cluster's element with the density it is tabulated at.
    """
    rows = [("outer", nucleus.name, nucleus.max_density_g_cm3) for nucleus in OUTER_CRUST]
    rows.append(("drip", "free neutrons", NEUTRON_DRIP_G_CM3))
    rows += [("inner", cell.name, cell.density_g_cm3) for cell in INNER_CRUST]
    return rows
