**Headline:** a body with more than a quarter of its surface under ice now shows an "Ice-capped" row on its sheet and page.

### Added
- GEN.1572: `hydrosphere.ice_capped` says a body is ice-capped when ice, sea ice and land ice together, covers more than 25 percent of its surface (`ICE_CAPPED_MIN_COVER`, registry row CC18, Marginal educated guess from physics, tweak level L0, a design edge Boss approved). The ice share is the larger of the cover the hydrosphere stored (`ice_cover_fraction`) and the share the body's own albedo law puts on a world of its mean surface temperature (`climate_branches.IceAlbedoLaw.ice_fraction`, the one GEN.1697 takes off the weathering land); on a liquid sea the edge falls at about 264 K. The flag is derived, not stored: the sheet (`describe/sheet.py`, `_ice_capped_provider`) adds an "Ice-capped" row to the water group of a planet or moon that is, and the pages and the sheet API show it. It changes no class letter and no liquid cover, and a body with no surface water (dry, vapor, a giant) never has it.

generation-output: unchanged (nothing stored, nothing generated differently: no column, no migration, no golden, no epoch, not a lock item)
generation-affects: none
