**Headline:** Habitability now rates a world by the climates people live in without machines, so a cold or hot place scores by how close it is to livable instead of by the clothing it needs.

### Changed
- **GEN.1677: the Temperature and UV edges of the PHI-6 score are condition-only.** The temperature edges are -120, -20, -10, 5, 31, 38, 122 and 200 C (Green from -10 to 38 C, Blue from 5 to 31 C) and the UV Green edge is 10 times Earth's index. Each too-cold or too-hot place is scored by its own band score, which replaces the flat Yellow credit. The gear thresholds stay as kit-list rows in `Bands.gear_temperature` and the protection limit table, apart from the rating.
- generation-output: changed
- generation-affects: fill
