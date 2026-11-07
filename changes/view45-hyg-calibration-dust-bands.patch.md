**Headline:** The sky calibration against the HYG catalogue now holds the default-dust sky to a band too, in the small and the Milky Way sized galaxy.

### Changed

- `src/tests/test_sky_hyg_calibration.py` (VIEW.45): each case also checks the share of the dust-free count that the sky through the galaxy's own gas-built dust keeps (`dusty_band`: 0.55 to 0.85 in the small galaxy, measured 0.693; 0.50 to 0.75 in the Milky Way sized one, measured 0.628), prints the default-dust tiers' reach against the dust-free tiers', and carries its own timeout so the suite's 900 seconds does not end a long run. The docstring has the new measurements.
- `reports/view45-hyg-calibration-dust-2026-11-07.md`: the measurements, the bands and what they say; one finding (a cosh overflow warning in the dust layer) left for its own item.

generation-output: unchanged
