**Headline:** The numeric precision audit says 64-bit floats are enough for the whole galaxy, and lists eight places to fix.

### Changed
- **Numeric precision audit (GEN.372).** `docs/design/numeric-precision-audit.md` gives, for every float the galaxy stores, computes or draws, how far from its origin it reaches, the spacing of a 64-bit float there (28 km at the Sun, 898 km at the 262 kpc reach, in parsecs), the smallest change that must survive, where the first error shows and a verdict. No stored value needs a wider type; the cases that need work (a local offset rebuilt from absolute coordinates, the field kernel on absolute metres, absolute 32-bit buffers in the browser, three phase formulas, the map tile root, a cell-key wrap at the reach edge, the region test) are listed with a named test each.
- **Two older documents corrected.** `orbital-updates.md` and `orbital-solvers-and-integrators.md` now say which unit the 32.8 km figure is for, and that positions are stored as DOUBLE milliparsecs, not integers.

TODO items resolved: GEN.372. Added: none.
