# planetgen/physics/nebula_columns.py

"""
The composition a nebula or a supernova remnant stores (GEN.1276, part g of GEN.1268, reports/gen1268-nebula-composition-2026-11-05.md).

Boss decided (8:49 PM EST on 2026-11-04) that the composition is **stored**, not derived at read time, so that as much as possible is searchable. This module turns
the stored row of a nebula (`nebulae`) or a remnant (`supernova_remnants`) into the values of the new columns, with the models of `nebula_composition` (the
gas, the dust, the phases), `planetary_nebula_ejecta` (a planetary nebula's progenitor and ejecta) and `remnant_composition` (a remnant's ejecta and swept-up
gas). It is the one place both the generator (right after a row is inserted) and the migration that fills the rows of an older galaxy call, so a re-run of
the backfill gives the same values: every draw is keyed on the object's ID (`row_key`).

**The key.** The object's `uid` as hex when it has one (a sector issues it when the row is inserted, and it is fixed by the seed and the place), else the
table and the name for an object saved on its own. An object that has a composition keeps it: a backfill fills only the rows that have none.

**The columns** (`COMMON_COLUMNS`, both tables; `REMNANT_COLUMNS`, remnants only), all nullable doubles, NULL for a row not filled yet:

* `gas_oh_dex`: [O/H] of the gas the object formed in (`nebula_composition.gas_oh`: the Milky Way gas gradient at its galactic radius plus the keyed scatter);
* `mass_fraction_h`, `_he`, `_z`: the hydrogen, helium and metal mass fractions of gas and dust together (add to 1); `mass_fraction_c`, `_n`, `_o`, `_ne`, `_s`,
  `_ar`, `_fe`: the mass fraction of each of the elements the pages show;
* `abundance_c` and the rest of the seven: 12 + log(X/H) of the gas phase (NULL where a hydrogen-free ejecta has none);
* `ionized_share`, `atomic_share`, `molecular_share`: the shares of the hydrogen mass in each phase (add to 1);
* `dust_to_gas`, `dust_to_metals`: dust mass over gas mass, and over metal mass;
* `grain_amorphous_olivine`, `grain_amorphous_pyroxene`, `grain_crystalline_silicate`, `grain_iron`, `grain_amorphous_carbon`, `grain_graphite`, `grain_pah`, `grain_sic`, `grain_ice`: the
  share of the dust mass in each grain species (add to 1 where there is dust; GEN.1376); the molecules of the ice mantle are `nebula_composition.ICE_MANTLE`, one table, not stored;
* `progenitor_mass_sol`: the initial mass of the star that made a planetary nebula (classes H to L) or a remnant; NULL for the other nebulae.

A remnant adds `ejecta_mass_sol` and `swept_mass_sol` (the mass thrown out by the explosion and the interstellar mass inside the blast wave, solar masses) and the
hydrogen, helium and metal mass fractions of each part (`ejecta_mass_fraction_h`, `_he`, `_z`, `swept_mass_fraction_h`, `_he`, `_z`), so the pages show the two as two
columns. The columns above hold the two together.

**The words.** `composition` (a nebula's plain-words sentence) and `dominant_species` are written from the numbers (`composition_sentence`, `species_summary`), not
copied from the family table: what the person reads is what is stored.

Satellite galaxies (GEN.1273): a nebula inside one takes the [Fe/H] of the satellite whose gas holds its place (`ComponentField.gas_fe_h`) and its gas [O/H] is drawn inside the Magellanic
range above it (`nebula_composition.gas_oh`, `SATELLITE_GAS_OFFSETS`, registry row NBC3), so a nebula of a Cloud is more metal-poor than one at the solar circle. A remnant has no stored
center yet and keeps the Milky Way gradient.
"""

import contextlib
import math

from planetgen import tuning
from planetgen.galaxy import scales, seed as galaxySeed
from planetgen.physics import nebula_composition as nc, planetary_nebula_ejecta as pne, remnant_composition as rc

ABUNDANCE_ELEMENTS = nc.ABUNDANCE_ELEMENTS
GRAIN_COLUMNS = {
    "amorphous_olivine": "grain_amorphous_olivine", "amorphous_pyroxene": "grain_amorphous_pyroxene", "crystalline_silicate": "grain_crystalline_silicate",
    "iron": "grain_iron", "amorphous_carbon": "grain_amorphous_carbon", "graphite": "grain_graphite", "PAH": "grain_pah", "SiC": "grain_sic", "ice": "grain_ice",
}
"""dict: The grain species of `nebula_composition.SPECIES` -> its column (GEN.1376 split the silicate and carbon columns of GEN.1276 and added SiC)."""

COMMON_COLUMNS = (
    "gas_oh_dex", "mass_fraction_h", "mass_fraction_he", "mass_fraction_z",
    *(f"mass_fraction_{e.lower()}" for e in ABUNDANCE_ELEMENTS),
    *(f"abundance_{e.lower()}" for e in ABUNDANCE_ELEMENTS),
    "ionized_share", "atomic_share", "molecular_share", "dust_to_gas", "dust_to_metals",
    *GRAIN_COLUMNS.values(), "progenitor_mass_sol",
)
"""tuple: The composition columns of `nebulae` and `supernova_remnants`."""

REMNANT_COLUMNS = (
    "ejecta_mass_sol", "swept_mass_sol",
    "ejecta_mass_fraction_h", "ejecta_mass_fraction_he", "ejecta_mass_fraction_z",
    "swept_mass_fraction_h", "swept_mass_fraction_he", "swept_mass_fraction_z",
)
"""tuple: The columns a remnant has beside `COMMON_COLUMNS`."""

COLUMNS = {"nebulae": COMMON_COLUMNS, "supernova_remnants": COMMON_COLUMNS + REMNANT_COLUMNS}
"""dict: Table -> its composition columns."""

UNITS = {
    **{column: "frac" for column in COMMON_COLUMNS + REMNANT_COLUMNS},
    **{f"abundance_{e.lower()}": "dex" for e in ABUNDANCE_ELEMENTS},
    "gas_oh_dex": "dex", "progenitor_mass_sol": "Msun", "ejecta_mass_sol": "Msun", "swept_mass_sol": "Msun",
}
"""dict: The unit each column is declared in (`db.column_units`)."""

PLANETARY = tuple(pne.CLASSES)
"""tuple: The classes with a progenitor drawn from the planetary nebula model."""

SPECIES_SHARE = 0.05
"""float: The share of the hydrogen a phase must hold to be named in `species_summary`. [Stand-in NBC14, tier T9, tweak L0: a wording cut, no physics.]"""

GRAIN_NAME_SHARE = 0.15
"""float: The share of the dust a grain kind must hold to be named in `species_summary`. [Stand-in NBC15, tier T9, tweak L0: a wording cut, no physics.]"""

_GRAIN_WORDS = {"amorphous_olivine": "amorphous olivine dust", "amorphous_pyroxene": "amorphous pyroxene dust", "crystalline_silicate": "crystalline silicate dust", "iron": "iron dust",
                "amorphous_carbon": "amorphous carbon dust", "graphite": "graphite dust", "PAH": "PAH molecules", "SiC": "silicon carbide dust", "ice": "ice mantles"}
_ELEMENT_WORDS = {"C": "carbon", "N": "nitrogen", "O": "oxygen", "Ne": "neon", "S": "sulfur", "Ar": "argon", "Fe": "iron"}


def in_galaxy(galaxy_seed, shape):
    """
    The context a composition is worked out in outside a sector's save: the galaxy's seed and thin-disk scale length bound, as `run_galaxy` binds them, because the gas gradient is
    quoted in the galaxy's own scale (`galaxy.scales`). Without it a plan's field rows and a migration's backfill would read the fallback galaxy's scale and differ from the numbers
    a sector's save stores. With no seed or shape (nothing planned) it changes nothing.
    """
    if galaxy_seed is None or shape is None:
        return contextlib.nullcontext()
    return galaxySeed.seeded(galaxy_seed, "nebula-composition", "context", disk_scale_length_pc=shape.disk_scale_length_pc, shape=shape)


FIELD_CLASSES = ("M", "N")
"""tuple: The classes of the molecular cloud field (`galaxy.nebula_field.TABLE_CLASSES`): their composition is keyed on the cloud's place, not on its object ID (GEN.1277)."""

FIELD_COLUMNS = tuple(column for column in COMMON_COLUMNS if column != "progenitor_mass_sol")
"""tuple: The composition columns of `nebula_field` (GEN.1277): `COMMON_COLUMNS` without the progenitor mass, which no molecular cloud has. FLOAT there, 4 bytes each."""


def field_key(center_pc):
    """The key of a field cloud (GEN.1277): its field cell and its center in whole milliparsecs, which the `nebula_field` row (at plan time, before any object ID) and the
    `nebulae` row the same cloud becomes when a sector stores it both know, so the two read one composition."""
    edge = tuning.NEBULA_FIELD_CELL_PC
    cell = "|".join(str(math.floor(c / edge)) for c in center_pc)
    return f"nebula-field|{cell}|" + "|".join(str(round(1000.0 * c)) for c in center_pc)


def row_key(table, row):
    """The key every draw of the object is keyed on: for a field cloud (class M or N, with a center) `field_key`, else its `uid` as hex, else the table and the name."""
    if table == "nebulae" and row.get("nebula_class") in FIELD_CLASSES and row.get("center_x_pc") is not None:
        return field_key((row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]))
    uid = row.get("uid")
    if uid:
        return bytes(uid).hex()
    return f"{table}:{row['name']}"


def _component_fe_h(row):
    """The [Fe/H] of the satellite the object lies in (GEN.1273), read from the components of the galaxy bound for the save (`galaxy.seed.current_galaxy_shape`), or `None` for an object
    of the Milky Way, an object with no center, or a galaxy with no satellites."""
    center = tuple(row.get(f"center_{axis}_pc") for axis in "xyz")
    components = getattr(galaxySeed.current_galaxy_shape(), "components", None)
    if components is None or None in center:
        return None
    return components.gas_fe_h(center)


def _gas_oh(table, row, key):
    radius = row.get("galactic_radius_pc")
    if radius is None:
        radius = scales.default_galactocentric_pc()
    return nc.gas_oh(float(radius), key, _component_fe_h(row))


def _common(made, oh, progenitor_mass):
    values = {"gas_oh_dex": oh, "progenitor_mass_sol": progenitor_mass}
    for column, value in made.values().items():
        values[column] = value
    for kind, column in GRAIN_COLUMNS.items():
        values[column] = made.species.get(kind, 0.0)
    return {column: (None if value != value else value) for column, value in values.items()}


def nebula_values(row):
    """
    The composition columns of a `nebulae` row (a dict with `nebula_class`, `density_cm3`, `radius_ly`, `galactic_radius_pc`, `hosts`, `uid`, `name`).

    Returns:
        dict: Column -> value, every column of `COMMON_COLUMNS`.
    """
    key = row_key("nebulae", row)
    nebula_class = row["nebula_class"]
    oh = _gas_oh("nebulae", row, key)
    ejecta, progenitor_mass, subtype = None, None, None
    if nebula_class in PLANETARY:
        found = pne.progenitor(nebula_class, nc.fe_h_for_oh(oh), key)
        ejecta, progenitor_mass = found.composition.elements, found.progenitor_mass_sol
    elif nebula_class == "X" and row.get("hosts"):
        subtype = str(row["hosts"]).split()[0]
    made = nc.makeup(nebula_class, oh, row["density_cm3"], row["radius_ly"], ejecta=ejecta, wolf_rayet_subtype=subtype)
    return _common(made, oh, progenitor_mass)


def remnant_values(row):
    """
    The composition columns of a `supernova_remnants` row (`remnant_class`, `progenitor_type`, `compact_remnant_kind`, `density_cm3`, `radius_ly`,
    `galactic_radius_pc`, `uid`, `name`).

    Returns:
        dict: Column -> value, every column of `COMMON_COLUMNS` and `REMNANT_COLUMNS`.
    """
    key = row_key("supernova_remnants", row)
    oh = _gas_oh("supernova_remnants", row, key)
    fe_h = nc.fe_h_for_oh(oh)
    found = rc.remnant_composition(key, row["progenitor_type"], row.get("compact_remnant_kind"), row["radius_ly"], row["density_cm3"], fe_h=fe_h)
    ejecta, swept, total = found.fractions()
    made = nc.makeup(row["remnant_class"], oh, row["density_cm3"], row["radius_ly"], ejecta=total)
    values = _common(made, oh, found.progenitor_mass_sol)
    values["ejecta_mass_sol"] = found.ejecta_mass_sol
    values["swept_mass_sol"] = found.swept_mass_sol
    for prefix, part in (("ejecta", ejecta), ("swept", swept)):
        x, y = part.get("H", 0.0), part.get("He", 0.0)
        values[f"{prefix}_mass_fraction_h"], values[f"{prefix}_mass_fraction_he"] = (x, y) if part else (None, None)
        values[f"{prefix}_mass_fraction_z"] = (1.0 - x - y) if part else None
    return values


def field_values(nebula_class, density_cm3, radius_ly, center_pc):
    """
    The `FIELD_COLUMNS` values of one field cloud (GEN.1277), as a tuple in that order: the same `nebula_values` a `nebulae` row of the cloud gets, from the cloud's own density and radius
    and its place (`galactic_radius_pc` is its distance from the galactic center, as `store._insert_field_nebulae` places it).
    """
    x, y, z = center_pc
    values = nebula_values({"nebula_class": nebula_class, "density_cm3": density_cm3, "radius_ly": radius_ly, "galactic_radius_pc": math.sqrt(x * x + y * y + z * z),
                            "center_x_pc": x, "center_y_pc": y, "center_z_pc": z, "name": "", "uid": None, "hosts": None})
    return tuple(values[column] for column in FIELD_COLUMNS)


def values_for(table, row):
    """The composition columns of a row of `table` (`nebulae` or `supernova_remnants`)."""
    if table == "nebulae":
        return nebula_values(row)
    if table == "supernova_remnants":
        return remnant_values(row)
    raise ValueError(f"no composition for table {table!r}")


def composition_sentence(values):
    """The plain-words composition of a nebula or a remnant, written from its stored numbers."""
    parts = [f"{100.0 * values['mass_fraction_h']:.1f}% hydrogen, {100.0 * values['mass_fraction_he']:.1f}% helium and "
             f"{100.0 * values['mass_fraction_z']:.2f}% heavier elements by mass"]
    phases = [(values["ionized_share"], "ionized"), (values["atomic_share"], "atomic"), (values["molecular_share"], "molecular")]
    shown = [f"{100.0 * share:.0f}% {word}" for share, word in phases if share >= SPECIES_SHARE]
    if values["mass_fraction_h"] > 0.0 and shown:
        parts.append("the hydrogen " + (shown[0] if len(shown) == 1 else ", ".join(shown[:-1]) + " and " + shown[-1]))
    parts.append(f"dust {values['dust_to_gas']:.2g} of the gas by mass")
    return "; ".join(parts) + "."


def species_summary(values):
    """What a nebula or a remnant is mostly made of, in a few words, written from its stored numbers: the hydrogen phases, helium, the three most
    plentiful of the seven elements the pages show, and the grains that hold a share of the dust."""
    names = []
    if values["mass_fraction_h"] > 0.0:
        for share, word in ((values["ionized_share"], "ionized hydrogen (H+)"), (values["atomic_share"], "atomic hydrogen (H I)"),
                            (values["molecular_share"], "molecular hydrogen (H2)")):
            if share >= SPECIES_SHARE:
                names.append(word)
    if values["mass_fraction_he"] >= 0.01:
        names.append("helium")
    elements = sorted(ABUNDANCE_ELEMENTS, key=lambda e: -values[f"mass_fraction_{e.lower()}"])[:3]
    names.extend(_ELEMENT_WORDS[e] for e in elements)
    names.extend(_GRAIN_WORDS[kind] for kind, column in GRAIN_COLUMNS.items() if values[column] >= GRAIN_NAME_SHARE)
    return ", ".join(names)


def sentences_for(values):
    """`(composition, dominant_species)` of a filled row, as stored."""
    return composition_sentence(values), species_summary(values)
