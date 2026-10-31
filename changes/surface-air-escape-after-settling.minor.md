**Headline:** An Earth-mass world in the habitable zone of an old Sun-like star keeps its air: the escape filter, which now sees the star, compares its loss to the pressure the air settled at.

### Fixed
- **The air fixed point lost the whole air to the escape filter (GEN.1199).** With the filter seeing the star, `surface_air.solve` handed it the pressure of the step it was on, which is the pressure of what stays, so the loss took its share of the air again at every step; and the reducing hydrogen mixtures of the first steps, which a settled air does not keep, took all of it. An Earth-mass world in the middle of the habitable zone of a 10 Gy Sun-like star ended airless in 3 of 4 seeds. The fixed point now settles without the escape, then compares the loss to that settled pressure (`mixture_draw.draw_mixture(escape_column_bar=...)`) and goes on with the loss in it; the same worlds keep an air in 8 of 8 seeds.
- Generator epoch bump (numbered at merge): the stored air, pressure, temperature, class and everything that reads them change for rocky bodies whose air loses gas.
generation-output: changed
