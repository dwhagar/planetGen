**Headline:** A tidal rheology module gives a body's Love number and quality factor from its mantle's temperature, melt and the tide's frequency, ready for the tide chain to use.

### Added
- **`physics/tidal_rheology.py` (GEN.954, stage 1a of the tide chain).** A pure function of mass, radius, tide frequency, mantle temperature, melt share and the body's own drawn mantle viscosity (Andrade or Maxwell rheology, Henning et al. 2009, Efroimsky 2012) that returns `k2` and `Q`. It gives the Moon its measured Love number and a Q of 61, Io `k2/Q` of a few hundredths in partial melt with no number fitted to Io, and the solid Earth a semi-diurnal Q of a few hundred. Nothing calls it yet; GEN.955 wires it into the tidal heat, the interior heat and the circularisation time.
- **Stand-in rows TR1 to TR9** in the gap-filling plan (section 6.23) and the stand-in list, each with its source tier and tweak level, and a design note, `docs/design/orbits/tidal-chain.md`.
