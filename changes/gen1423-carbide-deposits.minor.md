**Headline:** A carbon-rich rocky world (a star with C/O of 0.8 or more) now carries moissanite (silicon carbide) and cohenite (iron carbide) deposits; oxygen-rich worlds never do.

### Added
- **Two carbide deposit types (GEN.1423, `physics/carbide_deposits.py`).** After the ore pass each rocky planet, moon and rogue planet is assessed for `moissanite_crust` (the SiC of the carbon rock, grade 0.71 to 0.95, the share the bulk draw's graphite range leaves) and `cohenite_metal` (Fe3C in the carbon-saturated metal, grade 0.30 to 0.64). Both are rows of the existing deposit tables (no column, no migration); a body with silicate rock has both rows `absent_by_state`. Stand-ins CH1 to CH3 in section 6.228.

### Changed
- The stored deposits of every rocky body gain two rows, so the next epoch is a new one.

generation-output: changed
TODO items resolved: GEN.1423
