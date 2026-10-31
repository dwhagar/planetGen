**Headline:** Gas and ice giants are painted with belts, zones, wavy edges and storms, colored by their cloud chemistry.

### Added
- **Banded giants on the System Map, the 3D system view and the body page.** `static/giantpaint.js` replaces the old fixed stripe texture. The band
  count follows the spin and the size, the contrast follows the cloud-top temperature (ammonia-cloud giants show the most, methane ice giants are nearly
  plain), and the palette follows the stored methane, ammonia, water and hydrogen sulfide. Each giant looks the same on every visit (seeded from its ID).
  The look record gains a `gas` field. Display only: no schema change.
