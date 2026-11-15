**Headline:** A law for the light that reaches the water under ice and snow now exists (`physics/ice_light.py`), so the oxygen that life makes under a frozen lid can be counted by the light that gets through.

### Added
- GEN.1726: `ice_light.transmission` gives the share of the photosynthetic light that reaches the water under a lid, `(1 - albedo) exp(-k_ice h_ice) exp(-k_snow h_snow)`, with coefficients from measurements under Antarctic pack ice, Arctic first-year ice and a perennially ice-covered lake (sea ice k 1.26 per metre, clear lake ice 0.32, snow 30; bare sea ice albedo 0.6, snow 0.9). It reproduces Lake A (8.5 percent against 5.7 to 6.6 measured) and the under-0.1-percent of spring pack ice. `thickest_ice_m` and `deepest_snow_m` give the lid that still passes a given share of the light (about 6.6 m of sea ice for one part in ten thousand). Registered as stand-ins IL1 to IL6 in section 6.293 of the gap plan; design note `docs/design/habitability/light-through-ice.md`. Nothing reads it yet: the oxygen under ice is GEN.1723, and the thickness of a world's ice and snow is an open gap.

### Changed
- `scripts/build_references.py` no longer carries a person's address: the contact it sends to Crossref's lookups is read from `PLANETGEN_REFERENCES_CONTACT` and is empty by default. Five papers added to the reference list (R5264 to R5268).
