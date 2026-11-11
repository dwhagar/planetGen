**Headline:** A nebula's or remnant's page shows its composition at a glance, and a remnant sets what the explosion threw out beside the gas it swept up.

### Added

- The At a glance panel of a nebula and of a supernova remnant now lists the hydrogen, helium and heavy-element shares of its mass and its dust-to-gas ratio, read from the stored composition (GEN.1276), with their units. The full sheet below it already listed every stored value (GEN.1269 to GEN.1277); nothing is derived at read time, so no row carries a derived mark.
- A remnant's page has a two-column table, "What the explosion threw out and what it swept up": the mass and the hydrogen, helium and heavy-element shares of the ejecta beside those of the swept-up gas (`web/item_pages.py::ejecta_view`, the `remnant-ejecta` panel of `phenomenon.html`). It reuses the sheet rows, so a remnant that stores none of them has no table.
- Tests: `test_property_sheet.py` (the glance rows and the two columns from a stored nebula and remnant) and `test_web_system_phen.py` (the page).

generation-output: unchanged
generation-affects: none
