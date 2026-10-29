**Headline:** Generating a sector is about a third faster on a slow machine, and every stored value is exactly the same.

### Changed
- **Faster generation (PERF.201).** The gas equilibrium solver's matrix products, the Planck fraction in the infrared optical depth, the interior heat integration, the water continuum and band tables, and the band insolation of a tilt are now worked out with fewer, larger steps (or once instead of again). The stored values are bit for bit the same: tests compare the new sums with the old loops on random inputs, and the golden seed is unchanged.
