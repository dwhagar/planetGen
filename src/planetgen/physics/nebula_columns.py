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
* `grain_silicate`, `grain_iron`, `grain_carbon`, `grain_pah`, `grain_ice`: the share of the dust mass in each grain kind (add to 1 where there is dust);
* `progenitor_mass_sol`: the initial mass of the star that made a planetary nebula (classes H to L) or a remnant; NULL for the other nebulae.

A remnant adds `ejecta_mass_sol` and `swept_mass_sol` (the mass thrown out by the explosion and the interstellar mass inside the blast wave, solar masses) and the
hydrogen, helium and metal mass fractions of each part (`ejecta_mass_fraction_h`, `_he`, `_z`, `swept_mass_fraction_h`, `_he`, `_z`), so the pages show the two as two
columns. The columns above hold the two together.

**The words.** `composition` (a nebula's plain-words sentence) and `dominant_species` are written from the numbers (`composition_sentence`, `species_summary`), not
copied from the family table: what the person reads is what is stored.

Satellite galaxies: a nebula inside one has its gas read from the Milky Way gradient until GEN.367 places nebulae there (the model takes the component's [Fe/H]
and draws inside the Magellanic range, `nebula_composition.gas_oh`).
"""

from planetgen import tuning
from planetgen.galaxy import scales
from planetgen.physics import nebula_composition as nc
from planetgen.physics import planetary_nebula_ejecta as pne
from planetgen.physics import remnant_composition as rc

ABUNDANCE_ELEMENTS = nc.ABUNDANCE_ELEMENTS
GRAIN_COLUMNS = {"silicate": "grain_silicate", "iron": "grain_iron", "carbon": "grain_carbon", "PAH": "grain_pah", "ice": "grain_ice"}
"""dict: The grain kind of `nebula_composition.GRAINS` -> its column."""

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
    **{column: "ratio" for column in COMMON_COLUMNS + REMNANT_COLUMNS},
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

_GRAIN_WORDS = {"silicate": "silicate dust", "iron": "iron dust", "carbon": "carbon dust", "PAH": "PAH molecules", "ice": "ice mantles"}
_ELEMENT_WORDS = {"C": "carbon", "N": "nitrogen", "O": "oxygen", "Ne": "neon", "S": "sulfur", "Ar": "argon", "Fe": "iron"}


def row_key(table, row):
    """The key every draw of the object is keyed on: its `uid` as hex, else the table and the name."""
    uid = row.get("uid")
    if uid:
        return bytes(uid).hex()
    return f"{table}:{row['name']}"


def _gas_oh(table, row, key):
    radius = row.get("galactic_radius_pc")
    if radius is None:
        radius = scales.default_galactocentric_pc()
    return nc.gas_oh(float(radius), key)


def _common(made, oh, progenitor_mass):
    values = {"gas_oh_dex": oh, "progenitor_mass_sol": progenitor_mass}
    for column, value in made.values().items():
        values[column] = value
    for kind, column in GRAIN_COLUMNS.items():
        values[column] = made.grains.get(kind, 0.0)
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
