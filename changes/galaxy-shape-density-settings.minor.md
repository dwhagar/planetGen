### Added
- The galaxy shape settings (Generate page, New galaxy and Customize's Galaxy shape tab, and `planetgen plan`) now take the density on a spiral arm crest, the density between arms and a core density, with limits and the usual values (1.4, 0.6 and 0 for no core) filled in (ADM.49). The core is extra density at the very centre of the galaxy, on top of the bulge. The bulge amplitude is labelled Bulge density. The Galaxy Map's density model reads the new settings too, and an existing galaxy keeps its density.

### Changed
- `planetgen plan --arm-amplitude` is replaced by `--arm-density` and `--interarm-density` (the usual 1.4 and 0.6 are the old 0.4). Galaxy settings files saved with the old option name must be saved again.
