# Variable-star kinds (GEN.944, the cheap flags)

Code: `physics/variability.py`. Sheet: the derived "Variability" row of a star (`describe.sheet.VARIABILITY_KEY`), in the identity group. Tests: `tests/test_variability.py`, `tests/test_property_sheet.py`. Registry 6.138 (VR1 to VR3).

Nothing is stored and nothing in generation changes: the kind follows from the star's stored `state` (GEN.786), class, temperature and mass when the sheet is made, as the moon's longest eclipse is. A star that does not vary in a way its phase explains has no row.

| Star | Variability |
|---|---|
| state `cepheid` | classical Cepheid |
| state `asymptotic-giant`, cooler than 3,500 K | Mira or semiregular |
| class `PMS`, under 2 solar masses | T Tauri |
| class `PMS`, 2 to 8 solar masses | Herbig Ae/Be |

## Not here, and why

- **Carbon and S stars.** They have their own row, from the thermal-pulse count and the third dredge-up (GEN.1026, [carbon-stars.md](carbon-stars.md)).
- **Ap and Am stars.** They need a slow rotation and a magnetic field; the generator gives a star its spin and field at the end of its build, and has no chemically peculiar branch to put them on.
- **Semiregular variables of the red-giant branch.** The row marks the asymptotic-giant branch only.
