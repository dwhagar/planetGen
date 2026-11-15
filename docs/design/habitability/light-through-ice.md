# Light through ice and snow (GEN.1726)

Status: built (the law and its inverse); not yet read by anything. Code: `physics/ice_light.py`. Tests: `tests/test_ice_light.py`. Research: `reports/ice-light-and-lichen-cold-edge-2026-11-14.md`
(section 3) and `reports/life-adaptations-climate-2026-11-14.md` (GEN.1724). Boss, 9:17 PM EST 2026-11-14: "implement it, make it happen".

## 1. Why

A frozen-over ocean is not a dark one. Most of the Arctic Ocean's photosynthesis already happens in water under more than 85 percent ice cover, and perennially ice-covered lakes
keep oxygen-making mats going year round under 3 to 5 m of ice. The refuge term of GEN.1723 (the oxygen that life makes where light and any liquid water meet) therefore has to count lit
water under a lid, as the ocean's area times the light that gets through, and not as zero or as a pool. This module is the law for the light that gets through.

## 2. The law

`T = (1 - albedo) x exp(-k_ice x h_ice) x exp(-k_snow x h_snow)`

`T` is the share of the photosynthetically active light (400 to 700 nm) at the surface that reaches the water. It is the two-layer Beer-Lambert law of [R5264] (equation 3.1) with the
reflection at the top that [R5265] take off before they fit the same law to a lake.

| Row | Material | Value | Source | Tier |
|---|---|---|---|---|
| IL1 | k, sea ice (per m) | 1.26, range 0.8 to 1.5 | [R5264, R5267] | T1, L0 |
| IL2 | k, clear fresh lake ice (per m) | 0.32, down to 0.2 for thick perennial ice | [R5265] | T1, L0 |
| IL3 | k, snow (per m) | 30, range 10 to 100 (31.8 measured; 2.1 for one wind-packed layer; temperate 3 to 30) | [R5264, R5265, R5268] | T1, L0 |
| IL4 | albedo | bare sea ice 0.6 (0.5 to 0.7); lake ice 0.84; snow 0.9 (0.8 to 0.97) | [R5264, R5265, R5266] | T1, L0 |
| IL5 | albedo, melting ice | 0.4 (0.2 to 0.6) | [R5266] | T1, L0 |
| IL6 | snow cover depth | 0.03 m e-folding: how fast a thin snow takes over the albedo | no source; 10 cm of snow takes 95 percent of the light [R5268] | T7, L1 |

Snow is the main absorber, 30 times stronger per metre than sea ice and about 100 times stronger than clear lake ice, so the snow depth on a lid decides most of `T`.

## 3. Checks against measurements

| Case | Law | Measured |
|---|---|---|
| Lake A, 1.97 m of ice, snow removed [R5265] | 8.5 percent | 5.7 to 6.6 percent |
| Lake A, 2 m of ice and 0.5 m of wind-packed snow (k 2.1, albedo 0.97) | 0.55 percent | 0.45 percent |
| Antarctic pack ice in spring, 0.5 m ice and 0.2 m snow of albedo 0.8 [R5264] | 0.03 percent | under 0.1 percent |

The thin lid of a glacier melt lake (16 cm of ice) gives 33 percent against 21 to 23 percent in the one paper that measured it; that paper also says under 10 percent, a conflict the research could not settle.

## 4. How thick a lid still lets life work

The floor of GEN.1722 (0.04 umol photons m-2 s-1 as a daily mean, net growth) under an Earth-bright sky (a daily mean of about 400 umol m-2 s-1 in the lit season) is `T` of about 1e-4. The
thickest lid that keeps that: bare sea ice about 6.6 m; bare clear lake ice 23 m at k 0.32 and 30 m at k 0.25; snow alone about 0.23 m; 0.1 m of snow takes about 2.4 m of sea ice from the
depth. A world's own star and orbit set the surface light, so these depths move with them. They are examples only.

## 5. What it leaves out, and what it needs

* The coefficients are for the Sun's spectrum. Ice absorbs least at 400 to 500 nm and more toward the red, so a red star's light falls faster under a lid.
* The top 10 to 20 cm of sea ice attenuates more than the inside; it is folded into the albedo and the bulk k. Ice algae on the underside shade themselves (k 20 to 50); melt ponds let more through. Neither is counted.
* **The thickness of the ice and the snow is not built.** No law gives the generator the ice a world's ocean wears or the snow on it (a growth law for the ice from the temperature and the snowfall). Until one
  exists a caller must pass a thickness of its own. This is an open gap, filed with TODO.

## 6. Reuse

The generator's own ice albedos (`climate_branches.ICE_COVERED_ALBEDO`, `surface_air.ice_albedo_law`) are Bond albedos of a whole world, clouds included, so they are not the albedo of the ice surface to
photosynthetic light and are not reused; the stellar spectrum module is the caller's source of surface light (`biogenic_oxygen` takes it there already). Nothing else in the code holds an ice or snow light law.
