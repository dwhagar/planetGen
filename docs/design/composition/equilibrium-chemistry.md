# The equilibrium chemistry module (GEN.518)

Status: built 2026-10-24 by lane M (builder B10, the gas-phase part); merged into the one solver by GEN.672 (lane L, 2026-10-24). The solver is `planetgen/physics/gas_equilibrium.py` with the thermodynamics of `species_thermo.py` (`physics/data/nasa_glenn_species.json`); the regime edges of section 2 are `planetgen/physics/equilibrium_edges.py`, tested in `src/tests/test_equilibrium_edges.py`. The first module of GEN.518 (`equilibrium.py` with `thermo_nasa7.json`) held the same NASA data a second time and is deleted; the surviving solver reproduces the record's values (the five methane-carbon monoxide points within 1 K, the water and hydrogen edges within 0.1 percent), and the table below is the old module's.

## 1. What it is

A Gibbs-energy minimiser for an ideal-gas mixture at a temperature, a pressure and an element budget, written from the equations (the method of the NASA CEA program) and not copied from any code. The unknowns are the element potentials and the total moles, solved by Newton's method with a damped step; the first guess holds each element in the species that uses least of it, hydrogen first, and every case tried converges in five to six steps from 300 to 6000 K and 1 mbar to 100 bar. It reads no random number and sums with `math.fsum`, so the answer is the same on every machine. Gas only: no condensates (so no clouds), no ions (so no H-, no ionisation), ideal gas (about 100 bar and below), standard state 1 bar.

The thermodynamics are the NASA Glenn polynomials of `species_thermo` (99 species, 79 gases and 20 condensed phases; `DATA-LICENSE.md`); the solver takes the gases made of the elements given. The edges use the element budget of the independent computation (`equilibrium_edges.reference_budget`: hydrogen, helium, carbon, nitrogen, oxygen, sulfur and argon at the solar ratios of `mixture_giant.SOLAR_PER_HYDROGEN`, the metals scaled by the metallicity, carbon set to the C/O times the oxygen).

## 2. Test against an independent computation

The research record (`tunable-gaps-sources/group4-giants-bd-clouds.md` section 2.4) holds a table from another minimiser (easyCHEM's polynomials, scipy). This module agrees within 0.1 percent for the methane-carbon monoxide line and within 0.1 percent for the water and hydrogen lines (the recipe rows' anchor points were refreshed from it in GEN.723; they had been computed with the first module and sat up to 19 K from it):

| P (bar) | CH4 = CO (K), Z = 1 | water 90 percent (K) | water 50 percent (K) | H2 half atomic (K) |
|---|---|---|---|---|
| 0.001 | 726 | 2331 (2331) | 2567 (2567) | 2575 (2574) |
| 0.01 | 830 | 2542 (2542) | 2843 (2843) | 2885 (2884) |
| 0.1 | 966 (966) | 2793 (2792) | 3182 (3182) | 3279 (3278) |
| 1 | 1157 (1156) | 3096 (3095) | 3610 (3609) | 3795 (3794) |
| 10 | 1439 (1438) | 3471 | 4165 | 4501 |

(The record's value is in brackets.) At ten times solar metals the methane line falls by 40 to 270 K (892 K at 0.1 bar, 1052 K at 1 bar), as the record found. The tests pin these.

## 3. What changed in the recipes

Rows h1, h2 and h3 (the hot-Jupiter regimes) drew their edges from fixed temperatures (1500 and 2200 K). The record showed they are pressure-dependent outputs: at photospheric pressures the methane line sits near 970 K, not 1500 K, and 2200 K is only 1 to 4 percent water loss. The edges are now log-log interpolations of the computed lines in the pressure of the layer the regime is read at (the new stored quantity `photosphere_bar`, about 0.1 bar for the dayside emission): the methane-carbon monoxide line between h1 and h2, and the temperature at which a tenth of the water is destroyed between h2 and h3. The 100 K softness of the edges stays a design choice, named in the rows' limits. Metallicity is held at solar in the anchors; ten times solar moves the lower edge by up to 270 K at 100 bar and by 75 K at 0.1 bar, which the softness does not hide, so the dependence on metallicity is the next refinement (section 4).

## 4. What it does not do

- **Row n6 (the silicate-against-haze crossover).** The 950 K is the methane line plus the silicate sink depth; the module gives the methane line, but the sink depth needs a condensation model (clouds), which is not built. The row is unchanged.
- **Rows m2, m5, m6 (abiotic ozone, halogen gases).** The module has the species (HF, HCl, F, Cl) and can compute them; the inventories they need (the crust halogen budget, GEN.524) are not here.
- **Ions, condensates, non-ideal gas** (see section 1), and the metallicity and C/O dependence of the regime edges.
- **Who supplies `photosphere_bar`.** The recipes are not yet read by the generator, so the quantity is read from no body yet and no generated value changes; the wiring item draws or derives it.
