# The equilibrium chemistry module (GEN.518)

Status: built 2026-10-24 by lane M (builder B10, the gas-phase part). `planetgen/physics/equilibrium.py`, the data file `physics/data/thermo_nasa7.json` (`scripts/build_thermo_table.py`), tests in `src/tests/test_equilibrium.py`.

## 1. What it is

A Gibbs-energy minimiser for an ideal-gas mixture at a temperature, a pressure and an element budget, written from the equations of Gordon and McBride 1994 (the method of the NASA CEA program) and not copied from any code. The unknowns are the element potentials and the total moles, solved by Newton's method with a damped step; the first guess holds each element in the species that uses least of it, hydrogen first, and every case tried converges in five to six steps from 300 to 6000 K and 1 mbar to 100 bar. It reads no random number and sums with `math.fsum`, so the answer is the same on every machine. Gas only: no condensates (so no clouds), no ions (so no H-, no ionisation), ideal gas (about 100 bar and below), standard state 1 bar.

The thermodynamics are the NASA seven-coefficient polynomials of NASA TM-4513 (a United States government publication, so licence-compatible with the CC0 repository; `DATA-LICENSE.md`). Fifty-four species are kept: H, He, C, N, O, S, Ar, F, Cl, Na, K, Ti, V, Si and P compounds; the default solve uses the 39 of H, He, C, N, O, S and Ar. The element budget comes from a metallicity and a C/O on the Asplund 2021 scale of `tuning.STELLAR_SOLAR_LOG_EPS`, which is the reference scale and not a preset.

## 2. Test against an independent computation

The research record (`tunable-gaps-sources/group4-giants-bd-clouds.md` section 2.4) holds a table from another minimiser (easyCHEM's polynomials, scipy). This module agrees within 0.1 percent for the methane-carbon monoxide line and within 0.6 percent for the water and hydrogen lines:

| P (bar) | CH4 = CO (K), Z = 1 | water 90 percent (K) | water 50 percent (K) | H2 half atomic (K) |
|---|---|---|---|---|
| 0.001 | 725 | 2340 (2331) | 2571 (2567) | 2575 (2574) |
| 0.01 | 829 | 2553 (2542) | 2848 (2843) | 2885 (2884) |
| 0.1 | 966 (966) | 2806 (2792) | 3188 (3182) | 3279 (3278) |
| 1 | 1156 (1156) | 3111 (3095) | 3618 (3609) | 3795 (3794) |
| 10 | 1437 (1438) | 3490 | 4177 | 4500 |

(The record's value is in brackets.) At ten times solar metals the methane line falls by 40 to 270 K (892 K at 0.1 bar, 1052 K at 1 bar), as the record found. The tests pin these.

## 3. What changed in the recipes

Rows h1, h2 and h3 (the hot-Jupiter regimes) drew their edges from fixed temperatures (1500 and 2200 K). The record showed they are pressure-dependent outputs: at photospheric pressures the methane line sits near 970 K, not 1500 K, and 2200 K is only 1 to 4 percent water loss. The edges are now log-log interpolations of the computed lines in the pressure of the layer the regime is read at (the new stored quantity `photosphere_bar`, about 0.1 bar for the dayside emission): the methane-carbon monoxide line between h1 and h2, and the temperature at which a tenth of the water is destroyed between h2 and h3. The 100 K softness of the edges stays a design choice, named in the rows' limits. Metallicity is held at solar in the anchors; ten times solar moves the lower edge by up to 270 K at 100 bar and by 75 K at 0.1 bar, which the softness does not hide, so the dependence on metallicity is the next refinement (section 4).

## 4. What it does not do

- **Row n6 (the silicate-against-haze crossover).** The 950 K is the methane line plus the silicate sink depth; the module gives the methane line, but the sink depth needs a condensation model (clouds), which is not built. The row is unchanged.
- **Rows m2, m5, m6 (abiotic ozone, halogen gases).** The module has the species (HF, HCl, F, Cl) and can compute them; the inventories they need (the crust halogen budget, GEN.524) are not here.
- **Ions, condensates, non-ideal gas** (see section 1), and the metallicity and C/O dependence of the regime edges.
- **Who supplies `photosphere_bar`.** The recipes are not yet read by the generator, so the quantity is read from no body yet and no generated value changes; the wiring item draws or derives it.
