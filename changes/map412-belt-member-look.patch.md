**Headline:** A belt member now has a look record, so the pages and the 3D view can draw it: a lump or a ball in the colour of its taxon and makeup.

### Changed

- MAP.412: `bodylook.member_look(member)` gives the look record of a `belt_members` row: `kind` `belt_member`, a seed from its id, its radius and density, its taxon, and a colour (`member_grain_color`) that mixes the colour of the taxon's rock family (the same `ASTEROID_GRAIN_COLORS` the asteroid fields use; the dark taxa C, B, P and D carbonaceous, V basaltic, M metallic, the rest stony) with the metallic grey by its metal share and the icy pale blue-grey by its water and ice shares. `lumpy` is true when its roundness ratio at the median strength (`body_strength.ratio_for`) is under 1, so `lumpymesh.js` draws a lump; a body of Ceres' size is a ball. The info panel's account (`colour_account`) names the taxon and the shares.
- Reused: `seed_for`, `ASTEROID_GRAIN_COLORS`, `_srgb_byte`, `body_strength.ratio_for`. Not built here: the member page and the 3D view that call it (MAP.413, MAP.414), and members at or above the planet class floor, which draw through `planet_look` (GEN.1509 to GEN.1511). No query or stored value changes.

generation-affects: none
