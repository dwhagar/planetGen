# The air reads the system's late accretion (GEN.1493, design)

Status: design for review by TODO and Boss; nothing is built. Written by Lane 22 on 2026-11-09. Slice A (GEN.1453, `physics/influx.py`) is merged.

## The question

A rocky planet's air and ocean are drawn with one late-accretion number: `impact_history.late_accreted_fraction` gives 0.7 percent of the planet's mass times `mass^(-1/3)`, scattered
log-normally, and `surface_air._late_accretion` feeds it to `impact_history.veneer_balance` (what the impactors bring, less what they erode; GEN.739, GEN.819). That is one Solar-System
leftover for every system. A system with a heavy belt, a giant that sweeps the inner region clean, or a companion that stirs the reservoir has a different supply, and the air cannot know it,
because it is drawn when the planet is made and the belts and comets of the system come after.

## What the code does today (read, not assumed)

1. `Planet.__init__` calls `planetPhysics.calculate_atmospheric_conditions` (planet.py 468 and 500); `validation.py` calls it again where a planet is changed (183, 192, 605). This is the first and
   only draw of the late accretion (the gauss and the volatile scale come from the ambient stream and are recorded in `design`).
2. `StarSystem._add_belts` and `_generate_comets` run after the planet loop has settled (system.py ~420 and ~509). Habitable-world guarantees and `system.planets_regenerated` read the finished
   air, so the air must exist before the belts do.
3. `StarSystem._finish_planet` runs once the system is whole: `refresh_at_final_age` already remakes the hydrosphere, the ground and the score from the body's own stream (`draw.Stream` seeded
   from what is final about it), so a second pass that moves no draw of the run's stream is an established pattern. `Planet._surface_w` is kept (not stored) for it.
4. Slice A's `influx.planet_influx_kg_per_yr` is a present-day rate of dust and comets on a sheet. It is not the late accretion: that is the whole early supply, over the first hundreds of millions of
   years, when a belt of `M(t) = A / t` held several times its present mass. Belts store `primordial_mass_earth`, `cleared_mass_earth` and `scattered_mass_earth` (belt.py 126 to 183), so the early
   supply is on hand without a new stored value.

## Findings that shape the design

* **The slice A cross-section is the wrong geometry for a veneer.** `cross_section_share = R^2 / (4 a^2)` is the share of an isotropic stream that crosses the planet's orbit sphere once. For a
  planetesimal on a crossing orbit that is not lost, the chance of hitting the planet over gigayears is the sum over many encounters (an Opik/Wetherill collision probability with gravitational
  focusing, `R^2 (1 + v_esc^2/v^2) / (pi ...)` per orbit), several orders larger. By a back-of-envelope with the slice A law (a primordial reservoir of a few Earth masses, share 0.1, Earth's radius at 1 AU: 4.5e-10 of the stream) the direct share gives about 1e14 kg to an Earth-like planet against
  the published 0.2 to 0.7 percent of its mass (about 1e22 kg): some eight orders short. So the build needs the cumulative collision probability of a reservoir's planetesimals with the planet over the
  reservoir's decline, not the stream share. That is a new stand-in row set, and the calibration target is the published veneer range, not a Solar-System copy.
* **The giants cut both ways.** A giant between the reservoir and the planet ejects (slice A's shield, `scatter.THETA_EJECT`) but also scatters bodies inward (row IF6, not built). A first build keeps
  the shield only; the inward scattering is a named second step, because it is the part with no source in hand.
* **A lone star with the Solar System's architecture must stay inside the published range.** There is nothing fixed from the Solar System in the code (project rule), so the range is a test (a
  census over many Solar-System-like synthetic architectures), never an input.

## Options for the order

A. **Second pass at the end of the planet's making (recommended).** Keep the first solve with the present number as the provisional air the loop and the guarantees see. In `_finish_planet`, before
   `refresh_at_final_age`, work out the fraction from the system's reservoirs, and when it differs from the provisional one, solve the air again (`surface_air.solve`, same setting) with the new
   fraction and re-store it, then let `refresh_at_final_age` remake the hydrosphere, the ground and the score on it. The extra draws come from the body's own stream, so nothing of the run's stream
   moves. Cost: one more air solve per rocky planet (PERF: the solve is the dearest step; measure on a census, and skip the re-solve when the fraction moves under a stated tolerance). Downstream
   fields that read the first air before the pass (the class is re-applied in `_finish_planet`, so is the score; the habitable-world count in `habitable_request_status` is read after) follow.
B. **Make belts and comets before any air.** The air then reads the true reservoirs once. But the planet loop retries and guarantees on the finished air (habitable requests, forced classes), and
   belt birth reads the planet chain; the order would be circular, and every draw of the run's stream moves (a wide epoch, `--affects plan, scatter, fill`). Not recommended.
C. **Leave one number, make it depend on system properties.** A closed-form fraction from the star's metallicity, the giants' masses and the planet's place, with no reservoirs. Smaller, but it is the
   hand-fitted rule the item says to replace, and would not see a belt that is present or absent.

## Build in slices (each small, each its own PR)

0. **Calibration without wiring:** a pure `late_supply_kg(planet, reservoirs, giants, star)` (collision probability with focusing, the early reservoir from `primordial_mass_earth`, the slice A shield),
   a synthetic-architecture census against the published veneer range (0.2 to 0.7 percent for an Earth). Rows in GEN.725, the share and the collision constants honestly L3 or L4 (no source prints them).
   Lock-free: no output changes.
1. **Wire the second pass** (option A): the re-solve in `_finish_planet`, the first solve's gauss kept as the scatter about the new mean (so the distribution keeps its spread), a stored-field-free
   change. Lock item, `generation-affects: fill, population`; epoch with `--affects plan, scatter, fill` if the order of any draw of the run's stream moves, which option A is built not to do.
2. **Inward scattering by the giants and the companion's stirring** (rows IF6 and GEN.1492), after GEN.1492.
3. **Census and report:** the veneer by architecture (no giant, a Jupiter, a hot Jupiter, a companion), the share of planets whose air changed, and the cost of the extra solve.

## Done when

A Solar-System-like architecture (a census of synthetic ones, not the real one) gives a veneer inside the published range; a system with no belt and no comet gives a smaller one than a system with
a heavy belt and no giant; a giant between the belt and the planet lowers it by the shield; the cost of the second solve is measured and stated; stand-ins are registered with tier and tweak.

## Questions for Boss (plain words)

1. Is the published veneer range of the Earth (0.2 to 0.7 percent of its mass) the right target for the calibration, with the share and the collision constants honestly marked as fitted (L3 or L4)?
   Recommended: yes; nothing else measures it.
2. Is one extra air solve per rocky planet an acceptable cost for the end-of-making pass? Recommended: yes, with the skip when the fraction barely moves.

## Slice 0 result (built: `physics/late_supply.py`, rows LS1 to LS4 of section 6.250, census `tests/audit/late_supply_census.py`)

The law is the shared-sink form: each body takes the scattered part of a belt's early mass in proportion to `R^2 (1 + v_esc^2 / v_rel^2) / a^2`, so a giant next to a belt takes most of it and the planets
inside take what is left (the ejection is the denominator, not a second shield). With a stand-in inward share of 0.1 an Earth-like planet inside a belt of a few to tens of Earth masses reaches the
published 0.2 to 0.7 percent of its mass (test). In a census of 581 generated rocky planets of 0.3 to 3 Earth masses of single stars (`reports/gen1493-late-supply-census-2026-11-09.md`):

* the inward share at which the median planet that takes anything takes the mean of the present law is 0.105, so 0.1 stands (row LS1, L3);
* systems with a giant give a median 20 times lower than systems without one (the giant takes the share), a physical difference the single number never showed;
* 179 of the 581 planets (31 percent) take nothing, because no belt lies outside them: slice 1 needs an answer for these (a floor from planetesimals that never formed a belt, or a belt-free supply
  from the planets' own accretion tail); the single number gave them all something, so leaving them at zero would be a step backwards;
* the spread across planets (10th to 90th percentile a factor of about 1,500) is far wider than the log-normal scatter of the law (sigma 0.65 in ln, `impact_history.LATE_ACCRETED_LN_SIGMA`), which comes from belts that differ by orders of
  magnitude in mass; slice 1 either accepts it or narrows it with a floor, a call for Boss with the sample in hand.

## Slice 1 result (built: `generation/late_accretion.py`, `planets.remake_air`, rows LA1 and LA2 of section 6.251, census `tests/audit/late_accretion_census.py`)

- The pass runs once the system is whole, before the snow line pass, over the rocky planets of a star (not moons, giants or starless bodies). Each planet's percent is the supply of the system's belts and planets (`late_supply.veneer_fraction`), at least `FLOOR_SHARE` (0.05) of the law's mean for its mass, at most the old ceiling (3 percent).
- The air is made again by replaying `calculate_atmospheric_conditions` from the stream state saved at its last call (`planet._air_replay`, kept, not stored) in a private `draw.Stream`, so no draw of the run's stream moves. The air solves twice inside, each with its own draw of the late accretion; the replay holds both to the one percent, so a planet comes back as it was within about 1 percent in pressure and temperature when it drew nearly the same twice, and replaying twice with one percent is identical.
- Census (240 seeds, 246 planets of 0.3 to 3 Earth masses): 43 percent take the floor, 22 percent sit at the ceiling, the median is 0.13 of the law's mean (10th and 90th: 0.05 and 4.36). The pass is 18 percent of generation time (`reports/gen1493-late-accretion-census-2026-11-09.md`).
- Cost: the pass for every planet is 18 percent of generation time. Tried and dropped: a warm start of the air from the old solution (the fixed point has more than one answer for some planets, 3 of 10 trial planets landed elsewhere). Kept: a gate (`planets.LATE_ACCRETION_GATE`, row LA3): the air is left as it is when the new veneer moves the carbon and the nitrogen it left by less than half of what the volcanoes gave, and a planet with no air is not made again. The pass is then 8.8 percent; 201 of 246 planets are identical to the full replay, 13 differ by more than 1 percent (worst 99 percent in pressure, a climate branch). Gate 1.0 would be about 6 percent with a worst of 29 percent.
