# Aerosols: the titanium cold trap, row h4 (GEN.519, partly built)

Status: partly built 2026-10-23 by lane M (builder B11). Row h4 of the atmosphere recipes (`scripts/atmosphere_recipes_data.py`); tests in `src/tests/test_tuned_draws.py` and `src/tests/test_atmosphere_recipes.py`. Rows n6 and n7 are not built (section 2).

## 1. h4: TiO and VO on a hot Jupiter

The row was a function of the grain size of the cloud, `clamp((11 - (7.48 + 3.52 (log10 grain + 1) / 2)) / 4, 0, 1)`: a design formula that stood for the cold trap of titanium on the night side, whose depth depends on the vertical mixing (Kzz) that no population measures (two papers on one planet differ by 2.6 to 3.2 dex). Decision D5 of the gap plan drops the Kzz draw for a temperature threshold.

The row is now a step in the equilibrium temperature: TiO is present when `teq_k >= c_h4`, with `c_h4` drawn once per world uniform from 2350 to 2641 K. The bracket is the data: TiO is seen at 2641 K and not at 2210, 2250, 2260 and 2350 K (N = 5 planets; the threshold lies between the hottest non-detection and the detection). A uniform draw is the widest statement of an unknown threshold in a bracket, and a step is the form that adds no width the data do not give. Where the world has no equilibrium temperature the row does not apply.

Tests: no world has TiO at 2349 K, every world has it at 2641 K, and 25 and 50 percent of worlds have it a quarter and half way across the bracket (the uniform threshold).

## 2. What is not built

- **n6 (silicate against haze).** The row keeps `phi((teq - 950) / 100)`. The 950 K is the methane line plus the silicate sink depth: a model output that GEN.518 (the chemistry module) must compute before it can be tuned. The measured quantity (Fu 2017, N = 34: the water band amplitude rises 0.0008 +- 0.00016 per K) needs an amplitude model that links cloud cover to the band amplitude; it does not exist in the generator, so a target cannot be simulated.
- **n7 (haze yield).** The row stays a stand-in (`y` log-uniform 0.01 to 1). Every population statistic is degenerate with metallicity, Kzz and cloud cover, and the one published scatter ([R0516], about 0.5 scale height, N = 15) is in scale heights, not in a yield. It needs a haze amplitude model first.

Both are filed to the TODO thread as the remainder of GEN.519.

## 3. What it does not do

The recipes are not yet read by the generator, so no generated value changes and no epoch is taken.
