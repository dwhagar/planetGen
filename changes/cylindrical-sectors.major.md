### Changed
- **Galaxy sectors now sit on a cylindrical grid instead of spherical
  shells.** Each galaxy-placed sector is one cell of rings 11.5 ly wide
  around the galactic axis, layers 11.5 ly tall (layer 0 centered on the
  galactic plane) and wedge-shaped slots about 11.5 ly across, so sectors
  follow the flat disk instead of a ball. A sector's address is
  `(ring, layer, slot)`; systems and phenomena are placed inside the real
  cell, with local axes pointing outward, along the ring and north. See
  `docs/design/galaxy-coordinate-system.md`, "Cylindrical sector grid".
- `generate.py galaxy` takes `--ring I [--layer J] [--slot K]` in place of
  `--shell K [--slot N]`, and `--max-ring` in place of `--max-shell`.
  `generate.py plan` builds one band of layers per ring and takes
  `--max-ring`; `--workers` and `--chunk-size` are gone, since the build
  now takes about half a second. The Galaxy Map's copied commands use the
  new flags.
- The Galaxy page's 100 ly radial groups are now called Zones.
- Designations encode ring, layer and slot, so every sector gets a new one.

### Removed
- **Upgrading deletes every galaxy-placed sector.** Schema v31's migration
  deletes each placed sector together with its systems and phenomena, since
  shell addresses have no matching cell, and rebuilds the skeleton from the
  stored shape. Sectors that were never placed in the galaxy are kept.
  Regenerate the galaxy afterwards (for example `generate.py plan`, then
  `generate.py galaxy`), and take a backup first if you want the old data.
- The `sector_vertices` and `galaxy_shell_band` tables.
