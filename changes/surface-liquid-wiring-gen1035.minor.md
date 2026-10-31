**Headline:** Every rocky or icy world now says what liquid lies on its surface and how much of the surface it covers: water oceans, methane seas like Titan's, condensed seas, or the lava of a molten surface.

### Added
- **The surface liquid of a world (GEN.1035).** A new module (`physics/surface_liquid_world`) gives each planet and moon a stored `surface_liquid` (water, a hydrocarbon, a lava, and the other 50 liquids of the table), the share of the surface it covers and its mean depth. A water ocean keeps the cover and depth the hydrosphere found. A surface at or over the melting point of its ground's rock is a lava world (a magma ocean over the peridotite liquidus). Otherwise a cold body with an air gets the sea its ice supplies (a Titan-like body with 1 percent methane ice covers 1.7 percent of its surface, against Titan's 1.6) or the one its air condenses. Migration 0188 adds the three columns; the object page lists them under Water.

### Changed
- The stand-in registry has four new rows (LQ11 to LQ14), the share of an ice that lies on the surface, the lava rule, its pressure floor and the condensed sea mass.

generation-output: changed
