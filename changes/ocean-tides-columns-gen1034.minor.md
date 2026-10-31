**Headline:** A world with a sea now stores the height of its ocean tide, from its star and its moons, and the property sheet shows it.

### Changed
- **Ocean tides (GEN.1034, second part).** Every planet and moon with a surface ocean now has `tide_range_spring_m` and `tide_range_neap_m`, the open-ocean equilibrium tide in metres with every tide-raiser in line and with the largest against the rest, and `ocean_liquid` (water for now). An Earth with its Moon and Sun comes out near the published 0.5 to 0.8 m spring range, a sea with no moon shows only its star's tide, a world locked to its star shows none from it, and an ice-covered ocean shows 0 with no free surface. A coast can differ many times from the open-ocean value, because the basin's shape and resonance are not modelled; the property sheet says so beside the numbers. A new schema version adds the three columns; a world made before it keeps them empty until its system is generated again.

generation-output: changed
