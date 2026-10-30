**Headline:** The star density that rogue planets are drawn from can now be smoothed along the galactic orbit, so a rogue is no longer tied to the spiral arm or cluster of the star it left (not used by the generator yet).

### Changed
- **Smoothed rogue density (GEN.1115).** New `galaxy/rogue_density.py`: the density averaged round the galactic ring of a sector and over the radial and vertical spread an ejected planet gains (190 pc and 30 pc, from the drift module), with the few percent of recent ejections keeping the sector's own density. One function for the sector draw and the low-zoom map. Nothing calls it yet.

generation-output: unchanged
