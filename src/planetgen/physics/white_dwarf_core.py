# planetgen/physics/white_dwarf_core.py

"""
The core of a carbon-oxygen white dwarf and what a polluted one accreted (GEN.1370, N16 of the 8.1 real-material set)
=====================================================================================================================

Two things the white dwarf of the generator lacked:

1. **The core.** `white_dwarf_models.layer_mass_fractions` gives a carbon-oxygen white dwarf the layers of the downloaded sequences, which assume an
   even split of carbon and oxygen and no neon. A real core is oxygen-rich in the middle: the ten models of [R5238: Table 1] (progenitors of 1
   to 5 solar masses at Z = 0.01) have a central oxygen mass fraction of 0.62 to 0.79 that rises to a peak near 0.63 solar masses and falls
   to 0.62 at 0.88, and the two white dwarf pulsators fitted by asteroseismology have 0.61 and 0.84 [R5239]. `central_oxygen` reads the table
   (linear in the mass, held outside it) and holds the value to the asteroseismic range. The 22Ne is the CNO of the progenitor: each carbon,
   nitrogen and oxygen nucleus of its birth composition ends as 22Ne (`neon22_fraction`), 0.0152 for the proto-solar mix of the scale
   (`stellar_abundances.composition`; the 0.0137 of the photospheric solar abundances is the same sum without the 0.04 dex the proto-solar
   scale adds for the settling of the metals), and it scales with the progenitor's metallicity through that composition. `core_fractions` gives
   the three central mass fractions, the oxygen and carbon sharing what the neon leaves.
2. **The pollutant.** A polluted white dwarf shows the elements of the body it accreted [R1832]. `pollutant` takes the stored materials of a belt
   body (`small_body_composition.SmallBodyMakeup.components`, material id to mass fraction), reads each material's formula
   (`materials.json`; a few solid solutions and the alloys and organics have their own end members or shares in `tuning.py`), and returns the body's element mass
   fractions, the ratios to magnesium of the elements a white dwarf's spectrum shows (`OBSERVABLE`), and the class of the fragment: `core-rich`
   (at least `tuning.POLLUTANT_CORE_RICH_METAL_FRACTION` of the grains metal), `crust-rich` (a rock whose (Ca + Al) over Mg is
   `tuning.POLLUTANT_CRUST_RICH_FACTOR` times the proto-solar one), `volatile-rich` (water and ice at least
   `tuning.POLLUTANT_VOLATILE_RICH_FRACTION`) or `primitive` (the rest, which is most of them: more than 60 percent of the 202 systems of [R1832]).

Not here: the wiring into the white dwarf's stored layers and the pollution flag (a lock item that changes the output), the settling of
the elements out of the atmosphere (the spectrum shows the elements as they sink; [R1832] and the detection model `wd_pollution` treat it), and
the metallicity dependence of the central oxygen (the models of [R5238] are at Z = 0.01 only).
"""

from dataclasses import dataclass

from planetgen import tuning
from planetgen.physics import elements, materials
from planetgen.physics import small_body_composition as sbc
from planetgen.physics import stellar_abundances as sa

OBSERVABLE = ("O", "Mg", "Si", "Fe", "Ca", "Al", "Ti", "Cr", "Mn", "Ni", "Na", "K", "S", "C", "P")
"""tuple: The elements the spectra of the best-observed polluted white dwarfs give ([R1832], recalled set)."""

CORE_RICH, CRUST_RICH, VOLATILE_RICH, PRIMITIVE = "core-rich", "crust-rich", "volatile-rich", "primitive"
CLASSES = (PRIMITIVE, CORE_RICH, CRUST_RICH, VOLATILE_RICH)


class WhiteDwarfCoreError(ValueError):
    """A white dwarf whose core is not carbon-oxygen, or a composition with no carbon, nitrogen or oxygen."""


@dataclass(frozen=True)
class CoreFractions:
    """The central mass fractions of a carbon-oxygen core; they add up to 1."""

    oxygen: float
    carbon: float
    neon22: float


@dataclass(frozen=True)
class Pollutant:
    """
    What a polluted white dwarf would show of one accreted body.

    Attributes:
        kind (str): One of `CLASSES`.
        element_fractions (dict): Element symbol -> mass fraction of the body's grains, adding up to `1 - unassigned`.
        ratios_to_mg (dict): For each element of `OBSERVABLE` the body holds, its mass over the magnesium's; empty when the body holds no magnesium.
        unassigned (float): The mass fraction of materials with no formula, element share or end member (none for what a belt draws).
    """

    kind: str
    element_fractions: dict
    ratios_to_mg: dict
    unassigned: float


def central_oxygen(wd_mass_sol):
    """
    The central oxygen mass fraction of a carbon-oxygen white dwarf of `wd_mass_sol` solar masses: the table of `tuning.WD_CENTRAL_OXYGEN_BY_MASS_SOL`
    (linear in the mass, held outside its 0.5249 to 0.8779), held to `tuning.WD_CENTRAL_OXYGEN_RANGE`.

    Raises:
        WhiteDwarfCoreError: The mass is not one of a carbon-oxygen core (`stellar_abundances.white_dwarf_core`).
    """
    if sa.white_dwarf_core(wd_mass_sol) != "carbon-oxygen":
        raise WhiteDwarfCoreError(f"a {wd_mass_sol} solar mass white dwarf has a {sa.white_dwarf_core(wd_mass_sol)} core, not a carbon-oxygen one")
    points = tuning.WD_CENTRAL_OXYGEN_BY_MASS_SOL
    if wd_mass_sol <= points[0][0]:
        value = points[0][1]
    elif wd_mass_sol >= points[-1][0]:
        value = points[-1][1]
    else:
        value = points[-1][1]
        for (x0, y0), (x1, y1) in zip(points, points[1:]):
            if wd_mass_sol <= x1:
                value = y0 + (y1 - y0) * (wd_mass_sol - x0) / (x1 - x0)
                break
    low, high = tuning.WD_CENTRAL_OXYGEN_RANGE
    return min(max(value, low), high)


def neon22_fraction(birth):
    """
    The 22Ne mass fraction of a white dwarf core made from a progenitor of birth composition `birth` (a `stellar_abundances.StellarComposition`): the
    sum over carbon, nitrogen and oxygen of `tuning.WD_NEON22_PER_CNO` times the mass fraction times 22 over the atomic weight.

    Raises:
        WhiteDwarfCoreError: The composition holds no carbon, nitrogen or oxygen.
    """
    total = 0.0
    for symbol in ("C", "N", "O"):
        total += birth.elements.get(symbol, 0.0) * 22.0 / elements.element(symbol).atomic_weight
    if total <= 0.0:
        raise WhiteDwarfCoreError("the composition holds no carbon, nitrogen or oxygen to make 22Ne from")
    return tuning.WD_NEON22_PER_CNO * total


def core_fractions(wd_mass_sol, birth):
    """
    The central mass fractions of the carbon-oxygen core of a white dwarf of `wd_mass_sol` solar masses whose progenitor was born with `birth`: the 22Ne
    of `neon22_fraction`, and the rest shared by oxygen and carbon in the ratio `central_oxygen` gives (the tabulated fractions are of the
    carbon and oxygen alone).
    """
    neon = neon22_fraction(birth)
    oxygen = central_oxygen(wd_mass_sol)
    return CoreFractions(oxygen=oxygen * (1.0 - neon), carbon=(1.0 - oxygen) * (1.0 - neon), neon22=neon)


def _elements_of(material_id):
    """The element mass fractions of one materials-table row, or None when nothing states them."""
    shares = tuning.POLLUTANT_ELEMENT_SHARES.get(material_id)
    if shares is not None:
        return dict(shares)
    formula = tuning.POLLUTANT_FORMULAS.get(material_id) or materials.material(material_id).formula
    if not formula:
        return None
    atoms = elements.parse_formula(formula)
    mass = {symbol: count * elements.element(symbol).atomic_weight for symbol, count in atoms.items()}
    whole = sum(mass.values())
    return {symbol: value / whole for symbol, value in mass.items()}


def _crust_index(fractions):
    """(Ca + Al) over Mg by mass, or None without magnesium."""
    return (fractions.get("Ca", 0.0) + fractions.get("Al", 0.0)) / fractions["Mg"] if fractions.get("Mg") else None


def _proto_solar_crust_index():
    return (sa.PROTO_MASS_FRACTION["Ca"] + sa.PROTO_MASS_FRACTION["Al"]) / sa.PROTO_MASS_FRACTION["Mg"]


def accreted_taxon(members, u):
    """
    The taxon of the belt member a polluted white dwarf accreted: one of `members` (pairs of mass in kg and taxon letter, the stored members of the system's
    belts), chosen with a probability in proportion to its mass by `u`, a uniform number in [0, 1). `None` without members.
    """
    members = [(mass, taxon) for mass, taxon in members if mass > 0.0]
    if not members:
        return None
    target = u * sum(mass for mass, _taxon in members)
    running = 0.0
    for mass, taxon in members:
        running += mass
        if target < running:
            return taxon
    return members[-1][1]


def pollutant_of_taxon(taxon):
    """The `pollutant` of a body of belt taxon `taxon` (`small_body_composition.TAXON_DATA`): the taxon's own mix of the materials table."""
    return pollutant(sbc.TAXON_DATA[taxon].fractions)


def pollutant(components):
    """
    The elements a white dwarf would show of an accreted body whose grains are `components` (material id -> mass fraction, as
    `small_body_composition.SmallBodyMakeup.components`), and the class of the fragment (see the module).

    Raises:
        WhiteDwarfCoreError: No component has a known material or the fractions are empty.
        materials.UnknownMaterial: A component is not in the materials table.
    """
    total = sum(components.values())
    if total <= 0.0:
        raise WhiteDwarfCoreError("a body with no components has nothing to show")
    found, unassigned, metal, volatile = {}, 0.0, 0.0, 0.0
    for material_id, share in components.items():
        weight = share / total
        row = materials.material(material_id)
        group = sbc.group_of(row.id)
        if group == "metal":
            metal += weight
        elif group in ("water", "ice"):
            volatile += weight
        parts = _elements_of(row.id)
        if parts is None:
            unassigned += weight
            continue
        for symbol, fraction in parts.items():
            found[symbol] = found.get(symbol, 0.0) + weight * fraction
    if not found:
        raise WhiteDwarfCoreError("none of the components states its elements")
    magnesium = found.get("Mg", 0.0)
    ratios = {symbol: found[symbol] / magnesium for symbol in OBSERVABLE if symbol in found} if magnesium > 0.0 else {}
    if metal >= tuning.POLLUTANT_CORE_RICH_METAL_FRACTION:
        kind = CORE_RICH
    elif volatile >= tuning.POLLUTANT_VOLATILE_RICH_FRACTION:
        kind = VOLATILE_RICH
    else:
        index = _crust_index(found)
        kind = CRUST_RICH if index is not None and index >= tuning.POLLUTANT_CRUST_RICH_FACTOR * _proto_solar_crust_index() else PRIMITIVE
    return Pollutant(kind=kind, element_fractions=found, ratios_to_mg=ratios, unassigned=unassigned)
