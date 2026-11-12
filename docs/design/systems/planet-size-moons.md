# Moons of planet size around the giants (GEN.1157)

Boss's ask (10:26 PM EDT 2026-10-30): decide from physics how often a giant holds a moon the size of a planet and how large it can be, draw it with the same generator as every other body (no special case), and give it the same classes and life rules. Nothing from the Solar System is fixed content. Research: `reports/gen1157-planet-size-moons-2026-11-12.md`. Slice 0: `physics/planet_moons.py` (pure, unwired). Stand-ins: section 6.276 of the gap-filling plan (PSM1 to PSM4).

## 1. The physics in three routes

| Route | What it gives | Mass limit | Rate |
|---|---|---|---|
| Circumplanetary disk | regular moons, composition by the ice-line ramp | total of all its moons at most 6.4e-4 of the planet (usually 1e-4 to 2.5e-4), [R0642, R0643] | every giant with a disk; planet size needs a parent of several Jupiter masses ([R5229] restates the extension to ten) |
| Capture (binary exchange, or a planet caught during a migration instability) | one massive moon, any orbit, about half retrograde | a share of the planet; the candidate Kepler-1625 b-i is about 5.4e-3 [R5229]; the partner's mass sets the limit ([R5229: Eq. 1]) | rare: one candidate in 70 giants [R2365] |
| Giant impact | the rocky parents' moon | 0.006 of the parent and up | GEN.962, already built |

A disk cannot make an Earth-mass moon around a Jupiter: the whole disk holds about 0.2 Earth masses. Around ten Jupiter masses the ceiling is about 2 Earth masses, and the moons of that disk share it. That is why Solar System giants hold no planet-size moon while a rare super-Jovian or a captured partner can.

## 2. What the generator does today

`moon_composition` already makes the disk (`DiskSystem`, with the total from the same law), but the moon masses come from `planet_kind.draw_moon_mass_kg` with a flat ceiling of 1e-3 of a giant (`tuning.MOON_GIANT_PARENT_MAX_MASS_RATIO`). That is higher than the disk's top (6.4e-4), so a giant can already draw a moon above what its disk made; and there is no capture window above it. The wiring fixes both.

The Triton-class window of `irregular_population` (`exchange_chance`, GEN.1077) is a different thing: a planetesimal binary gives an ice or rock body of 1e22 kg or so, and it holds no gas giant. The planet-size capture here starts at 0.1 Earth mass (6e23 kg), so the two never give the same moon; the wiring checks that a planet-size capture is not also drawn by that window.

## 3. Proposed wiring (changes output, rides a train)

1. In `generate_moons`, cap a disk-origin moon at `planet_moons.disk_moon_ceiling_kg` (and the disk's sum at the same figure), replacing the flat 1e-3 for giants.
2. After the disk moons, one `planet_moons.draw_exchange_moon(planet.mass, draw)` per giant; a hit adds a moon of that mass with origin `exchange` (a new value of `moon_composition.ORIGINS`), placed anywhere inside `moon_orbit_bounds_km`, with the captured inclination law (isotropic, about half retrograde, [R5229]), and made by the ordinary `Planet` generator, so its class, air and life come from its own conditions.
3. Epoch with `--affects`; the census shows how many giants gain a planet-size moon (about 1.4 percent by design) and how many lose a heavy disk moon.

The hosting gate (`physics/hosting_gate.py`) still decides which giants can hold moons at all, so a hot giant keeps none.

## 4. What is not modelled

The binary encounter itself (the bound of Eq. 1 of [R5229] needs the ejected partner), tidal circularization of the captured orbit (the orbit here is drawn, not integrated), and the survival of a moon through the planet's migration. Each is a named stand-in or a later slice.
