**Headline:** A sector fill runs as stages that each draw from their own random stream, and a scattered pulsar gets its planets only when its sector is filled.

### Added
- **The fill's stages (GEN.251).** `generate_sector` runs the stages of the phenomena catalogue: A the pre-placed stars, B the scattered rows, C the system pass, D the cell's other systems, E the per-star kinds (one stream per kind), F the cloud field and the nebulae grown around hot stars, G the neutron stars and black holes below the scatter's mass cut, H the fast-star flags and the nucleus. Each draws from `draw.Stream(f"{base}:{stage}")`, `base` being one number drawn from the sector's own stream when the fill begins (`planetgen.generation.fill_stages`), so adding a kind to a stage, or a stage, moves no draw of another stage. The sector keeps what each stage made by kind (`sector.stage_counts`) and the debug log prints it for every sector.
- **The system pass (stage C).** A scattered or below-cut neutron star or black hole rolls for planets once, on `draw.Stream(f"{seed}:system-pass")` (the scatter row's seed; a below-cut remnant's: its sector's below-cut seed and its place in the draw), so the answer is the same whenever and in whatever order the sector is filled. A millisecond pulsar holds planets 0.7% of the time and becomes a star system at its own point (`SpaceSector.replace_phenomenon_with_system`); no other compact remnant holds any. Nothing exists around a scatter row before its sector is filled. The stellar companions of those objects stay in the scatter (GEN.262) until the fill builds them from the stored members (GEN.264).
- **`tuning.SYSTEM_PASS_ODDS`** holds the chances of what an object can hold beyond its stars, by kind (a neutron star's by pulsar type, a black hole's by mass class): the millisecond pulsar's planets (0.7%) and the rogue planet's moon (20%). The phenomena registry names the model (`remnant-planets`) for neutron stars.

### Changed
- **The rogue planet's moon chance moved from `ROGUE_PLANET_MOON_CHANCE` to `SYSTEM_PASS_ODDS["rogue-planet"]["moons"]`** (the same 20%).
- **A seeded sector fills differently from before.** Every stage now draws from its own stream, so the same galaxy seed gives different systems and phenomena than the single stream did.

TODO items resolved: GEN.251. Added: none.
