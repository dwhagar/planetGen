**Headline:** Nebulae on the Galaxy Map fade out at the edge by family: a reflection nebula over the longest way, a supernova shell crisp.

### Changed
- **MAP.142 (the Galaxy Map's part).** The strength of a nebula's points now follows `smoothstep(FORM_LEVEL - w / 2, FORM_LEVEL + w / 2, density)` with a width `w` for each family (`galaxy/nebula_volume.EDGE_WIDTH`: reflection 0.60, diffuse 0.50, emission 0.40, dark 0.30, planetary 0.15, supernova remnant 0.10 of the density's span, from the "Appearance by class" table), so a nebula has no hard outline and its extent is the half-strength contour, the density the model calls inside a nebula, whatever its width. A planetary nebula and a remnant keep their shell profile, so the families differ by more than color: shells are rings, solid nebulae fade over their own width. The widths are display settings, not measurements. The Sector Map's and the nebula page's meshes keep their edges for now.

generation-output: unchanged

TODO items resolved: none (MAP.142 stays open for the Sector Map and the nebula view).
