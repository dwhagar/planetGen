**Headline:** The dust layer no longer makes numpy warn of an overflow for stars far above the galactic plane; its values are unchanged.

### Changed

- `physics/sky_photometry.cosh_squared` (VIEW.48): `cosh(x) ** 2` as numpy gives it where it is finite and `inf` beyond the float range (|x| above about 355.6), without the `RuntimeWarning: overflow encountered in cosh` and `in square`. The dust layer `1 / cosh(z / h) ** 2` of the gas-built profile (`dust_model.DustProfile.av_per_pc`) and of the fixed Milky Way law (`sky_photometry.dust_av_per_pc`) use it, so far above the plane the layer is exactly 0 as before. The result is identical for every input, nan and infinity included (a test compares the two over a grid and the edge values).
- `src/tests/test_dust_model.py`: the layer far above the plane in both laws under warnings as errors, and `cosh_squared` against numpy.

generation-output: unchanged
