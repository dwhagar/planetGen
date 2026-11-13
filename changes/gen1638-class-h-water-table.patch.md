**Headline:** A desert world with a little surface liquid can now be class H, which no world could be before.

### Fixed
- **GEN.1638: class H was unreachable.** H needs groundwater, and the aquifer model only stores groundwater for a world with no liquid at its surface, so a world with a brine seep or a small ocean never had any. `class_domains.derived` now counts a water table under the land of a `surface ocean` or `brine seep` world with an ocean cover over 0 and under 0.2 (H's own cut, held equal by a test). The stored `groundwater` column and every habitability score that reads it are left as they were, so a world with liquid at its surface still draws its light and takes its surface dose. In 160 census systems no class moved; H still needs the oxygen of GEN.1581 and GEN.1584.
generation-output: changed
generation-affects: fill
