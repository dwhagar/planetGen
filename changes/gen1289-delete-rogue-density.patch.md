**Headline:** The unused smoothed star density for rogue planets is deleted.

### Changed
- **GEN.1289: `galaxy/rogue_density.py` is deleted.** It was written for a virtual donor that drew a sector's rogues from a smoothed star density. A rogue now comes from a stored donor that threw it out, placed where it lands, and the top-up makes up what a sector lacks, so nothing reads the density. It goes with its tests, stand-in section 6.102 (rows RSM1 to RSM4) and its register lines. Generation output is unchanged.
