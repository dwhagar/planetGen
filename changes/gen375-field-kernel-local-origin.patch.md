**Headline:** Paths and pulls in a sector far from the galaxy's centre are as exact as in one beside it.

### Fixed
- The gravity kernel (`physics/field_kernel.py`) measures from a local origin: the start of a path, or the first point of a pull. A step is no longer added to a position 262 kpc from the centre, where the numbers are 1,049 km apart, so the low bits of every step were lost and the loss built up along the path. A path at the reach now matches the same path beside the centre to the spacing of the far numbers; before, it was off by about ten spacings.
- Stored paths are left out of the fingerprints and the golden files, so nothing in them moves. A sector's paths come out the same except for the last digits.

generation-output: unchanged
TODO items resolved: GEN.375
