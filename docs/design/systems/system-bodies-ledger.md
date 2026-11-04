# The system body ledger and the barycentre that includes the belts (GEN.1066)

Status: pure module built (`physics/system_bodies.py`, lock-free, nothing calls it yet); the view over the typed tables is the wiring step. Written 2026-11-04 as the lane reported. Follows belt-evolution.md section 5 (Boss, 4:19 PM EDT 2026-10-30) and unified-system-barycentre.md.

## 1. The question and the answer

Planets, moons, belts, comets and rings each live in their own typed table, so no query answers "what does this system weigh" or "where is its centre of mass with the belts in it". The decision asked for: **a view, not a table.** One generated list of rows, built from the typed tables when asked, with no new storage and no second copy to keep in step. A stored table would duplicate every mass and every orbit and would go stale the moment one of its sources changed; the rows are cheap to build (a system has at most a few hundred stored bodies, because a belt stores at most 300 L1 members), and the view is the single place the sums are made.

## 2. The row

`id`, `parent`, `kind`, `mass_kg`, `radius_km`, `density_kg_m3`, `a_au`, `e`, `i_deg`, `status`, `tier`, and `position_au` (relative to the parent's centre, as the typed tables store it).

| kind | comes from | parent | notes |
|---|---|---|---|
| star | `stars` | none (root) or a star | a companion star is a body about the primary here; close pairs keep the tree convention in `system_tree.py` |
| planet | `planets` | star | status `planet` or `dwarf planet` (belt-evolution.md 3.2) |
| moon | `moons` | planet, or a belt member that has a moon | |
| belt | `asteroid_belts` | star | a container: mass zero, its mass is the sum of its children |
| belt_member | the L1 rows of a belt (GEN.1069) | belt | tier `L1`, status `minor body` or `dwarf planet` |
| extended | the belt's remainder (L2 and L3, `belt_ledger.Ledger.remainder_kg`) | belt | tier `L3`; one point at the belt's centre |
| comet | `comets` | star | |
| ring | the ring system of a planet | planet | one extended mass at the planet's centre |

The view reads exact column names when it is written (the typed tables are the source of truth; this note does not repeat them). `check` makes the list a tree: one root, parents that exist and are of an allowed kind, no cycle, finite non-negative masses, and no mass on a belt itself.

## 3. The sums

* **Mass is conserved** (belt-evolution.md section 3). `subtree_mass_kg(belt)` is the L1 members plus the extended remainder, and `belt_mass_residual(rows, belt, belt_mass)` is that against the belt's stored mass, so a ledger that does not add up is a number. `belt_ledger` already guarantees at least five percent of the belt stays in the remainder.
* **Positions are stored relative, derived absolute** (unified-system-barycentre.md section 2). A body's place in the system frame is the sum of the offsets down its chain of parents (`absolute_positions_au`), so a moon is its planet's place plus its own offset and nothing above needs to know.
* **The barycentre is `sum(m r) / sum(m)`** over every row with a mass, belts and rings included (`barycentre_au`); the central object sits at minus that (`central_offset_au`). For one planet it is the two-body formula, for a moon the moon's share of the total mass. This is the "reflex" the planets columns store today for planets alone, now with the belts in it.
* **The belt test.** `belt_wobble_au` is the distance between where the barycentre puts the star with the belts and without them. A round, centred belt gives zero; one heavy member on one side gives `m a / (M + m)`. The belts' wobble of a generated star is therefore a measurable number, and the Solar System's (the Sun's wobble from the belts is a test target, never an input) is the first check once the view is wired.

## 4. What is not here

* The view over the typed tables, and where it is read (the system page, the map, the N-body export): the wiring item, which needs the belt member rows (GEN.1069) and the L2 on-demand draw to exist.
* A time series: positions are one instant; a caller loop advances the orbits (`relative_position_au` turns elements and a mean anomaly into a position).
* Close binaries as the origin of two stars, and moons of moons: the first keeps its own convention, the second is excluded by the no-satellite line of the unified rule (a belt member that has a moon is the only exception, as small-body binaries exist).
* Per-body gravitational interaction: the ledger is for sums, not integration.

## 5. Tests

`tests/test_system_bodies.py`: the tree checks and their reasons, chained positions, subtree and total mass, the barycentre against the two-body formula and with a moon, a round belt against a heavy member (the wobble), the belt residual, and elements to positions.
