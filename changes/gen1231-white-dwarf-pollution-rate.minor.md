**Headline:** A polluted white dwarf now shows the rate at which rock falls on it: its belt's own when the star's mass loss threw enough of the belt inward, else the typical rate of a polluted white dwarf.

### Changed
- `physics/belt_history.rate_on_white_dwarf` (GEN.1231) gives a polluted white dwarf a rate: `pollution_rate_g_s` of the mass its belts' star threw inward when that is enough to pollute one (`is_polluted`), else `TYPICAL_POLLUTION_RATE_G_S` (1e8 g/s, stand-in 6.117 BHS6, tier T2, tweak L0). Which stars are polluted does not change: the flag stays the drawn one (25 to 50 percent of the cool white dwarfs), because the belts of the white dwarf systems built so far threw nothing inward (none of 300 had both a belt and a giant). The star's sheet gains a "Rock falling on the star" row.
- `asteroid_belts` gains `thrown_inward_kg` (DOUBLE, NULL; schema v265, migration 0265): the mass the star's mass loss threw inward from the belt (`remnant_orbits`), NULL for a belt of no remnant system and for a row made before. Nothing is backfilled; a regenerated galaxy fills it. No other stored value changes.

generation-affects: none
