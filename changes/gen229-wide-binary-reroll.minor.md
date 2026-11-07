**Headline:** `planetgen bring-to-epoch` now finds the wide binaries that break the stability rule and makes only those systems again from their own seeds, saying which it redid.

### Added
- **Re-roll of unstable wide binaries (GEN.229).** A wide (S-type) binary whose planets or stars' outer bodies break the stability rule GEN.228 put on new ones is drawn again as a wide binary for its sector, from a seed of its own (the galaxy seed and the system's Object ID), keeping its row, name, Object ID and place (`galaxy/wide_binary_reroll.py`). The console names each system, the bodies that broke the rule and what the re-roll made. It is a stage of `planetgen bring-to-epoch` (`--no-wide-binaries` skips it, `--dry-run` only lists). A system that hosts facilities, is built around a placed bright star, or is in an edited sector is left alone and named.
- `generation/sector_context.py`: the stellar context of a sector, shared by the API's system routes and the re-roll.

### Fixed
- **A regenerated system kept the old orbit tree.** `store.replace_star_system_content` left `orbit_nodes` and `dark_members` of the old system in place, so a regenerated wide binary loaded with the old tree's stability limits. The tree now goes with the new stars.

generation-output: unchanged

TODO items resolved: GEN.229.
