# Nebula composition: what every nebula, remnant and cloud is made of

**Status (2026-11-05):** built as the model of GEN.1269 (part a of GEN.1268); the writes are GEN.1276 (nebulae and remnants) and GEN.1277 (field clouds), the planetary-nebula ejecta GEN.1270, the remnants' ejecta GEN.1271. The research and the arguments are in `reports/gen1268-nebula-composition-2026-11-05.md` (sections 2, 3, 4 and 11); this note holds what the code does and where it differs from that report.

**Boss's decisions** (8:49 to 8:55 PM EST on 2026-11-04): the gas has its own metallicity gradient, a separate stand-in; the Magellanic gas takes the widest range of the two published estimates; the dust survivals of report section 11; a Wolf-Rayet ring and a planetary nebula are enriched in helium and nitrogen; each field cloud has its own varied composition keyed on its cell; the composition is stored for nebulae, remnants and field clouds (stars stay derived); the pages show everything.

## 1. The module

`physics/nebula_composition.py` stores nothing. `makeup(nebula_class, gas_oh_dex, density_cm3, radius_ly, ejecta=None, wolf_rayet_subtype=None)` returns a `NebulaMakeup`:

* **Composition** of gas and dust together (`x`, `y`, `z`, and the mass fraction of each element), from the thin-disc birth composition (`stellar_abundances.composition`) whose [O/H] is the place's, so a nebula and a star born in it share one chemistry.
* **Dust**: the dust-to-gas ratio of [R5102] (`dust_to_gas_law`), times the survival of the class and held to what the metals allow; the dust by element; the gas left over; the share of the dust in silicate, iron, carbon, PAH and ice.
* **Phases**: the ionized, atomic and molecular shares of the hydrogen (a table by class; for the cold classes F, M, N, P and Q the molecular share of the neutral gas follows the column density and metallicity, [R5104]).
* **Abundances** of the gas phase, 12 + log(X/H), of C, N, O, Ne, S, Ar and Fe (`nan` for a hydrogen-free ejecta).
* `values()` is the dict the composition item stores.

`gas_oh(galactic_radius_pc, key, component_fe_h=None)` is the place: the gas gradient (0.034 dex per kpc, zero at the reference radius) with 0.14 dex of scatter from the nebula's key, or, in a satellite component, its [Fe/H] plus an offset drawn inside the published ranges of the two Clouds (report section 4). `field_cloud_density_cm3(class, key)` gives a field cloud its own density inside the class range.

## 2. Where this differs from the report

* The grains are five kinds (magnesium silicate, **iron** as its own grain, carbon, PAH, ice). The report's silicate carried its iron; with a fixed iron content a metal-poor gas ran out of iron first, so the iron is a grain of its own and the dust follows what the gas has. A solar diffuse nebula's dust is 0.83 of the law, the rest not being in the metals; the law value is kept in `dust_to_gas_law`.
* The caps (`MAX_IN_DUST`) keep every element's dust under its total, so the dust of a metal-poor nebula is whatever its metals allow.
* Class S uses the sourced shocked-remnant survivals (silicate 0.39, iron 0.60); H II regions and planetary nebulae keep 0.9 of silicate and iron grains and 0.5 of the carbon half. The fresh-ejecta dust of R, T, U and W is a placeholder 0.2 until GEN.1271 writes the ejecta rule of report section 11.
* A Wolf-Rayet ring mixes the host's surface (`stellar_composition.wolf_rayet_surface`) into the swept-up gas at `WIND_SHARE_OF_RING` (0.3, argued).
* The stand-ins are rows NBC1 to NBC13 of section 6.179 of the stand-in registry.

## 3. Open

Krumholz, McKee and Tumlinson's form is checked only in its chi coefficient and its half-molecular column (34 against the 27 of the paper's summary); [R5105] coefficients, the depletion of the Clouds, and the dust of old remnants stay to be read (GEN.1274). The progenitor of planetary nebulae and the ejecta of remnants enter through the `ejecta` argument (GEN.1270, GEN.1271); without it those classes read as interstellar gas and say so in `open_items`.

## 4. What is stored (GEN.1276, GEN.1277)

**Nebulae and remnants** (GEN.1276): `physics/nebula_columns.py` turns a row into the columns of `nebulae` and `supernova_remnants` (DOUBLE, searchable through the API object lists), at insert and in migration 0212. `composition` and `dominant_species` are written from the stored numbers.

**Field clouds** (GEN.1277, design note). About 163,000 rows in `nebula_field` (classes M and N) carry only a cell, a class, a radius and a place, about 70 bytes a row. Two ways to add the composition were weighed:

* *Stored float columns* (chosen): the 28 composition columns of `nebulae` without the progenitor mass, each FLOAT (4 bytes, 7 digits stored, read back as 6 significant digits by the text protocol, enough for fractions down to 1e-7 and for [O/H] to 1e-6 dex). 112 bytes a row, about 18 MB for 163,000 rows. Each is a plain column a range search or a sort reads, with a unit in `column_units`, and none needs decoding.
* *Packed integers* (one 16-bit word per value, or all of them in a BINARY blob): 56 bytes a row, about 9 MB. The saving is 9 MB, and the price is that nothing could be searched or sorted in SQL without a generated column per value (which brings the bytes back), and the 16 bits resolve a fraction only to 1.5e-5, which is 15 percent of a sulfur or argon fraction (1e-4). Rejected.

Table growth: from about 11 MB of data to about 29 MB, and a plan writes it in the same pass that already draws the cloud (no second pass).

The key of a field cloud is its cell and place (`nebula_columns.field_key`: the cell index and the centre in milliparsecs), not its object ID, because the row exists at plan time and has none until a sector stores it. The `nebulae` row the cloud becomes is keyed the same way (a class M or N nebula with a centre, `row_key`), so the two read one composition; the stored nebula keeps full doubles and the field row the FLOAT of the same number. The cloud's density is the one its `Nebula` was drawn with (`cell_draw`), which differs from cloud to cloud inside the class range, so gas, dust and phases differ with it and with the place's gas. `field_cloud_density_cm3` is kept for a cloud the draw does not give back (a table planned under another seed), in migration 0215 only. No progenitor mass is stored: no class M or N cloud has a progenitor.

