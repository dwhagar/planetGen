**Headline:** A comet on a closed orbit never leaves its star: its eccentricity is cut so the aphelion stays inside the star's Hill sphere.

### Changed
- **An elliptical comet's eccentricity is drawn from its period class's range cut at the value that puts the aphelion at its star's limit** (`generation/comet.py`, `generation/system.py`). The long-period class reaches e = 0.999 and, for a light star in the inner galaxy, an aphelion past the Hill sphere (a 0.08 solar mass star: 6159 AU against 5815 AU); GEN.182 says a closed orbit never leaves its star. `Comet` takes a `ceiling_au`: the star's Hill sphere against the galaxy for a single star or a close pair (`validation.orbit_ceiling_au`), the tighter stable limit for a star of a wide pair, which the system's re-roll already enforced. A range that lies wholly above the cut gives the cut itself. A parabolic comet is unbound on purpose and reads no ceiling. Draws that the cut does not reach are unchanged.
- The generator epoch is 37.
generation-output: changed

### Added
- `tests/test_comet_data.py`: the bound test hands the comet its Hill sphere (300 draws for each of three masses), a ceiling cuts the range and leaves a wide one alone, a parabolic comet reads none. `tests/test_systems.py`: every closed comet orbit of a generated system is inside its star's orbit ceiling.

TODO items resolved: GEN.591. Added: none.
