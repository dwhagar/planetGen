# The air reads the system's late accretion (GEN.1493, design)

Status: slices 0 and 1 are built and on main (see "What is on main" below); slices 2 and 3 are designed in the last section and not built. Written by Lane 22 on 2026-11-09; the last section by Lane 06 on 2026-11-18. Slice A (GEN.1453, `physics/influx.py`) is merged.

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

## Slice 0 result (built: `physics/late_supply.py`, rows LS1 to LS4 of section 6.254, census `tests/audit/late_supply_census.py`)

The law is the shared-sink form: each body takes the scattered part of a belt's early mass in proportion to `R^2 (1 + v_esc^2 / v_rel^2) / a^2`, so a giant next to a belt takes most of it and the planets
inside take what is left (the ejection is the denominator, not a second shield). With a stand-in inward share of 0.1 an Earth-like planet inside a belt of a few to tens of Earth masses reaches the
published 0.2 to 0.7 percent of its mass (test). In a census of 581 generated rocky planets of 0.3 to 3 Earth masses of single stars (`reports/gen1493-late-supply-census-2026-11-09.md`):

* the inward share at which the median planet that takes anything takes the mean of the present law is 0.105, so 0.1 stands (row LS1, L3);
* systems with a giant give a median 20 times lower than systems without one (the giant takes the share), a physical difference the single number never showed;
* 179 of the 581 planets (31 percent) take nothing, because no belt lies outside them: slice 1 needs an answer for these (a floor from planetesimals that never formed a belt, or a belt-free supply
  from the planets' own accretion tail); the single number gave them all something, so leaving them at zero would be a step backwards;
* the spread across planets (10th to 90th percentile a factor of about 1,500) is far wider than the log-normal scatter of the law (sigma 0.65 in ln, `impact_history.LATE_ACCRETED_LN_SIGMA`), which comes from belts that differ by orders of
  magnitude in mass; slice 1 either accepts it or narrows it with a floor, a call for Boss with the sample in hand.

## Slice 1 result (built: `generation/late_accretion.py`, `planets.remake_air`, rows LA1 to LA3 of section 6.255, census `tests/audit/late_accretion_census.py`)

- The pass runs once the system is whole, before the snow line pass, over the rocky planets of a star (not moons, giants or starless bodies). Each planet's percent is the supply of the system's belts and planets (`late_supply.veneer_fraction`), at least `FLOOR_SHARE` (0.05) of the law's mean for its mass, at most the old ceiling (3 percent).
- The air is made again by replaying `calculate_atmospheric_conditions` from the stream state saved at its last call (`planet._air_replay`, kept, not stored) in a private `draw.Stream`, so no draw of the run's stream moves. The air solves twice inside, each with its own draw of the late accretion; the replay holds both to the one percent, so a planet comes back as it was within about 1 percent in pressure and temperature when it drew nearly the same twice, and replaying twice with one percent is identical.
- Census (240 seeds, 246 planets of 0.3 to 3 Earth masses): 43 percent take the floor, 22 percent sit at the ceiling, the median is 0.13 of the law's mean (10th and 90th: 0.05 and 4.36). The pass is 18 percent of generation time (`reports/gen1493-late-accretion-census-2026-11-09.md`).
- Cost: the pass for every planet is 18 percent of generation time. Tried and dropped: a warm start of the air from the old solution (the fixed point has more than one answer for some planets, 3 of 10 trial planets landed elsewhere). Kept: a gate (`planets.LATE_ACCRETION_GATE`, row LA3): the air is left as it is when the new veneer moves the carbon and the nitrogen it left by less than half of what the volcanoes gave, and a planet with no air is not made again. The pass is then 8.8 percent; 201 of 246 planets are identical to the full replay, 13 differ by more than 1 percent (worst 99 percent in pressure, a climate branch). Gate 1.0 would be about 6 percent with a worst of 29 percent.

## What is on main (read 2026-11-18)

| Piece | Where | State |
|---|---|---|
| Early supply as a shared sink (rows LS1 to LS4) | `physics/late_supply.py` | built; every planet and giant is a sink, so a giant's ejection is the denominator |
| End-of-making replay with the system's number (rows LA1 to LA3) | `generation/late_accretion.py` (`apply`, `give_first_air`), `planets.remake_air`, `build.air_replay` | built; private stream, no draw of the run's stream moves; gate `LATE_ACCRETION_GATE` |
| Present-day dust and comets on a sheet (rows IF1 to IF10) | `physics/influx.py`, `generation/cosmic_dust.py`, `describe/influx_row.py` | built; display only, not the veneer |
| Giants' swept disk and companion stirring of belts | `physics/scatter.py` (`REACH_ZONE_WIDTHS` 3, GEN.1788), `physics/belt_forcing.py` (GEN.1492), `physics/belt_origin.py` (`stirred`, BG8) | built for the belt's mass and orbit; **not yet read by the veneer** |

So the veneer today is `INWARD_SHARE` (LS1, a constant 0.1, T6, L3) times each belt's **primordial** mass, shared among the bodies. Two things the system does to a reservoir never reach it:

1. A giant beyond the belt does more than take a share. Its chaotic zone throws a part of the belt's mass onto crossing orbits (`scatter_shares`: the bound, scattered part) or out of the system (the ejected part). Today the belt's `scattered_mass_earth` and `cleared_mass_earth` are stored and drawn on the belt page, and the veneer ignores them.
2. A companion star or a close giant stirs the belt (`stirred`, `belt_forcing.combined` raises its eccentricity). A stirred belt grinds faster and sends more inward; a quiet one sends less. LS1 is the same 0.1 for both.

Comets are the other gap. `cosmic_dust` gives them a present-day dust rate, but `late_accretion.percent_for` builds its reservoirs from belts only, so an outer disk that made comets and no belt gives its planets only the floor.

## Slices 2 and 3 (design, not built)

### Slice 2: the inward share reads what the system did to the reservoir

**Change.** Replace the constant in LS1 by a per-reservoir share, keeping the single row (the Boss tweak rule of 2026-11-16: modify the existing stand-in, no modifier layered on top):

    inward_share(belt) = LS1_BASE + LS1_GAIN x scattered_share(belt) x (1 - ejected_share(belt)) + LS1_STIR x stirred(belt)

- `scattered_share` is `belt.scattered_mass_earth / belt.primordial_mass_earth`, the part the giants threw onto bound crossing orbits; `ejected_share` is already inside `scatter_shares` (`ejected / swept`), so a heavy giant (large Safronov number) sends most of what it scatters out of the system instead of in. That is the physical difference between a Jupiter and a Neptune: the lighter giant keeps what it scatters.
- `stirred` is 0 or 1 today (BG8). The companion's stirring (`belt_forcing.combined`) already moves the belt's eccentricity, so a continuous form is possible: the stirring term goes as the belt's eccentricity over the unstirred width eccentricity. First build uses the stored flag, because it is already stored and a test can pin it.
- `LS1_BASE` stays 0.1, so a quiet belt behind no giant gives what the census calibrated in slice 0 (median 0.105 of the law's mean). The two gains are new constants, registered next to LS1 as LS1a and LS1b, tier T6 (physics-shaped, no source prints the number), tweak L4 (fitted, no source), each with a range test so no gain can lift a Solar-System-like synthetic architecture out of 0.2 to 0.7 percent of an Earth's mass.
- Comets: add the comet reservoirs of the system to `percent_for` (`influx.Reservoir("comet", ...)`, mass as the comet bodies' summed mass, ring at the perihelion), scaled by the same `inward_share` of the disk that made them. A system whose outer disk made comets and no belt then stops being a floor-only system.

**What it costs.** No new air solve. The number only feeds the existing replay, so the cost stays the 8.8 percent the gate left (slice 1 result). The share is evaluated per planet per belt, a sum over a handful of rows.

**Lock and epoch.** The veneer moves for every system with a giant or a stirred belt, so the air, the water and the class of some rocky planets move: a lock item, `generation-affects: fill, population`. No draw of the run's stream moves (the replay is private), so no new epoch is needed beyond the merge-train golden. Pass `--affects plan, scatter, fill` only if the comet reservoirs are read from a stage that runs before the planet loop, which this design avoids by reading them in `late_accretion.apply`, after the system is whole.

**Tests.** Pure-function tests on the share (no giant gives base; a Neptune keeps more of what it scatters than a Jupiter; a stirred belt gives more than a quiet one), the Solar-System-like range test over the synthetic census, and the wiring test that the replay reads the new share.

### Slice 3: census and report by architecture

`tests/audit/late_accretion_census.py` grows an architecture column (no giant, a Jupiter, a Neptune, a hot Jupiter, a companion star, a comet-only disk) and the report `reports/gen1493-veneer-by-architecture-<date>.md` gives, for each, the median veneer, the share of planets at the floor and at the ceiling (today 43 and 22 percent), the share of planets whose air changed, and the generation time with the pass. The done-when of this note then closes.

### Questions for Boss (plain words)

1. A belt behind a giant that throws things inward should feed the planets more, and a belt behind a very heavy giant should feed them less, because that giant throws most of what it scatters out of the system. No source prints how much more or less, so two new numbers would be tuned so a Solar-System-like system stays in 0.2 to 0.7 percent of an Earth's mass, and marked as fitted (L4). Is that acceptable? Recommended: yes.
2. Should comets count as a source of the veneer, the same way belts do? Recommended: yes; today a system that made only comets leaves its planets at the floor.
3. Today 43 percent of rocky planets sit at the floor and 22 percent at the ceiling, a wider spread than the single number gave. Slice 3 will report the spread by architecture; should the floor and ceiling stay as they are until that report is in? Recommended: yes, change nothing until the census is in front of you.

## Slice 3 result (baseline before slice 2; `tests/audit/late_accretion_architecture_census.py`, report `reports/gen1493-veneer-by-architecture-2026-11-18.md`)

In 240 Sun-like systems 71 percent of rocky planets sit on the floor (51 percent for default stars), 65 points of it because no belt lies outside the planet; every planet of a system with a giant (14) is at the floor. The median veneer is 0.04 percent of the planet's mass against the published 0.2 to 0.7 percent. So the design of slice 2 above changes: the reservoir to read is the one the system **had** (the belt's primordial mass, what the giants cleared and scattered, the outer disk that made comets), not only the belts that remain, or a giant that cleared its belt leaves the planets inside it at the floor. If the primordial mass of a fully cleared belt is not stored, storing it is the first step of slice 2. No comet-only system appears in 480 seeds (GEN.1788 rates), so the comet source needs a forced fixture.
