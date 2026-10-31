**Headline:** A small body is now an Asteroid or a round Class D or C world by whether its own weight beats its strength, not by how much ice it holds.

### Changed
- Class D is now an airless, round body under 0.1 Earth masses, rock or ice; Class C is an airless body at 0.1 Earth masses or more. The 1 percent water rule is gone, and moons use the same classes (GEN.1160).
- A body whose roundness ratio is under 1 is a new object type, the Asteroid. It has no class letter (its class is empty in the data, its pages say Asteroid), and it can be a moon. The ratio comes from the roundness engine (`physics/roundness.py`) at the body's own drawn strength (GEN.1159).
- The round mass of a dwarf planet candidate, the round test of the body status words and Class D of a belt member now all read the roundness engine: an icy body turns round at a few 1e19 kg and rock at about 2e20 kg at the median strength. The typed round masses (BR9, DPM1) and the compaction-taper roundness (BR24) are retired; a Kuiper-type disc now lists a few hundred dwarf candidates at the median strength, about 20 at the strongest.
- The class sweep and the temperate-ocean table count Asteroids as their own row.

generation-output: changed
