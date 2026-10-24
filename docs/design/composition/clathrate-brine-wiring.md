# Clathrate and brine limits in the generator (GEN.762)

Status: built 2026-10-24 (lane "Source and tune atmosphere draws"). Code `src/planetgen/physics/clathrate_brine.py` (`hydrate_hold_pa`, `draw_brine`, `water_activity`), `physics/atmosphere.py` (rule 5), `physics/hydrosphere.py` (the ice-covered branch); tests `src/tests/test_clathrate_wiring.py`. The curves and eutectics themselves are GEN.471 ([body-composition-research.md](body-composition-research.md)).

## 1. The atmosphere: cages in the ice skin

An icy body (water at least `ICY_MIN_WATER_FRACTION`, 1 percent of its mass) whose surface is under 273.15 K holds the guest gases CO2, CO, N2, CH4 and H2S in the water cages of its ice. A guest above its dissociation pressure `ln P = A0 + A1/T` forms hydrate until its pressure falls to the curve, so its demand is the excess over the curve. A pure cap at the curve would remove Titan's nitrogen (the N2 curve is 108 Pa at 94 K), so the cages are a finite store: one guest per six waters in a skin of ice whose depth is drawn per body, log-uniform 1 m to 1 km. The guests share the cages, so a small store takes the same share of every guest's excess. The mole count of a partial pressure is `p / (M g)`, `M` the mean molar mass of the stored gases (a helium or ammonia balance is not counted).

The held gas leaves the air as the condensate does in rule 4: the stored partial pressures and the total surface pressure are lower. Nothing new is stored (no column, no migration); the trapped amount is not kept.

## 2. The hydrosphere: brine under the lid

A body whose surface water is under `BRINE_MAX_WATER_FRACTION` (2e-5, the line of the chloride brine class) has its salts concentrated. When its surface is under the ice point, the base of the lid melts at the eutectic of a drawn perchlorate (NaClO4, Mg(ClO4)2 or Ca(ClO4)2, equal weights, 198 to 236 K, with the module's spread) instead of 273.15 K, less the pressure term as before. A thinner lid, or none, follows, and a brine warmer than its eutectic is a `surface ocean`, not covered. The ocean's water activity is that salt's figure where the source gave one (Mg 0.555, Ca 0.515; Na has none, so it still draws from the class range).

## 3. Stand-ins and limits

- The skin depth range (1 m to 1 km) has no measured figure: a stand-in, bracket as stated, replace when a regolith or kinetics source is read.
- Equal salt weights: a stand-in; Mars-like oxidised bodies favour perchlorate, and chloride eutectics (the class is named for them) are not in the data file yet.
- The brine applies only below the chloride-brine water line; Europa-like oceans with sea-salt freezing-point depression of a few K are unchanged.
- NH3 monohydrate is in the curve file but not a stored species, so it is not held.
- It changes generated values (the atmosphere and the ocean of cold icy and cold dry bodies), so it takes a generator epoch.
