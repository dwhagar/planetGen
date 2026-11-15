**Headline:** The generator can now say how thick the ice and the snow on a world's water are, so the light that gets through a frozen ocean or lake can be worked out for a real world.

### Added
- GEN.1728: `physics/ice_growth.py` gives the ice and the snow on a world's water from its mean surface temperature, the swing of its year, the starlight it absorbs and the heat the water gives up: the ice a freezing season grows (Stefan's law with the snow and the air in the way, the freezing degree days of the year in closed form), the thickness of a lid that never thaws (a heat balance against the water's flux and the wear of the surface, 2.9 m in the Arctic standard case against the 2.88 m of the model it comes from), and the snow of the year (the Earth's share of the starlight that goes into the water cycle, about two thirds of it kept as on the Arctic ice). Stand-ins IG1 to IG8 in section 6.298 of the gap plan; design note `docs/design/habitability/ice-and-snow-growth.md`.

### Changed
- `continental_weathering.latent_j_m3` is the latent heat of its energy limit, factored out of `runoff_m_yr` with the same value. Nothing reads the new law yet (the refuge term of GEN.1723 will): no body, no stored value and no page changes.
