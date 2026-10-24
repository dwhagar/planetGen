# Spherical sectors: meridians for the block ladder

Research only, 2026-10-15. Boss's goal (06:46Z): as many common meridians as possible, as few disjointed changes in meridian as possible. Builds on [spherical-sector-coordinates.md](spherical-sector-coordinates.md) (the plan, with wedge counts in multiples of 6). No code changes. Boss decides the questions in section 5.

Evidence tags: [C] computed with the scripts in `/mnt/project-files/research/scripts/spherical-sectors/` (`06_` to `08_`, results in `results-ladder-*.txt`), on the real default density model and outline. Cells are the outline's cells, shells 4 and out.

## 1. Answer in brief

1. **The plan's rule gives the best cell shape and the most meridian changes.** Wedge counts in multiples of 6 keep arcs within 0.1% of 4 pc but change the count in nearly every band: 157 changes per shell on average inside the populated window (95th percentile 770), and the equator count differs in every shell, 3,763 changes along the plane, none of them nested. Only the six pole meridians are common to everything [C].
2. **A nested chain removes every broken meridian at a price in arc accuracy.** Counts that each divide the next (6, 18, 36, ... up to the equator count) never end or shift a meridian: a meridian starts at one latitude and runs to its mirror in the south. Measured: 0 non-nested changes, 1.8 changes per shell inside the populated window (95th percentile 8), 2,745 of 3,761 shells (73%) need just one count over their whole populated window, and only 51 changes of equator count along the plane. The price: arcs 0.73 to 1.27 of 4 pc for the middle 98% of cells, 89% of cells within 10% of 4 pc, 97% within 25%, worst 0.58 to 1.73, and the cell count rises 2.1% (24.43 billion against 23.93) [C].
3. **Where the nested error falls.** It is two things. Along the equator, the allowed counts are `6 * 2^a * 3^b`; there are 34 of them between 630 and 24,000 wedges and neighbours differ by up to 18.5%, so an equator arc is off by up to 9.3% [C]. Toward the poles, each halving step lets the arc swing 0.71 to 1.41 [C]. In the thin disk of a far shell the populated window is only a few degrees wide, so the whole window uses one count and the error is the equator's alone.
4. **Radial sharing comes with it.** With `6 * 2^a * 3^b` counts, a run of shells shares one equator count, so all their meridians are common; the 34 counts mean about 50 changes along the plane. At a change, neighbouring counts keep a median 12.5% of their meridians (every eighth), not none [C].
5. **There is a middle choice, not measured**: fine multiples of 6 inside the populated disk window and a nested chain only above it (one non-nested change per shell, at the window edge).

## 2. Terms

- **Common meridian:** a wedge boundary present in every band of a stretch of the shell (or of a run of shells). More is better for drawing block edges and for picking blocks that do not cut a cell.
- **Disjointed change:** a step in wedge count where an old meridian ends or moves, so the new count is not a multiple of the old (`W(b+1)` not a multiple of `W(b)`). A nested step only starts new meridians.
- **Populated window:** the bands of a shell that the outline keeps. It is 805 of 1,570 northern bands at shell 1,000 and 184 of 4,712 at shell 3,000.

## 3. What was measured

### 3.1 Inside one shell (latitude) [C]

Per northern hemisphere, pole cap to equator, `06_latitude_chains.py`:

| Shell | Rule | Equator count | Changes in populated window | Not nested | Meridians common to the whole populated window |
|---|---|---|---|---|---|
| 100 | multiples of 6 | 630 | 103 | 103 | 6 |
| 100 | nested chain | 648 | 5 | 0 | 6 |
| 1,000 | multiples of 6 | 6,282 | 333 | 333 | 6 |
| 1,000 | nested chain | 6,144 | 6 | 0 | 12 |
| 3,000 | multiples of 6 | 18,852 | 6 | 6 | 6 |
| 3,000 | nested chain | 18,432 | 0 | 0 | 18,432 (all) |

At shell 3,000 the populated window is 184 bands wide and the chain gives it a single count: all 18,432 meridians are common, arcs 1.02 of 4 pc. This is what the cylinder's aligned columns gave for free.

### 3.2 The whole galaxy [C]

`07_wholegalaxy_wedge_structures.py`, outline cells, shells 4 and out:

| | Multiples of 6 (plan) | Nested chain, `6 * 2^a * 3^b` equator |
|---|---|---|
| Cells | 23.933 billion | 24.433 billion (+2.1%) |
| Arc / 4 pc, 1st / 50th / 99th percentile | 0.999 / 1.000 / 1.001 | 0.73 / 0.99 / 1.27 |
| Cells within 10% / 25% of 4 pc | 100% / 100% | 89% / 97% |
| Changes per shell in populated window, mean (95th) | 157 (770) | 1.8 (8) |
| Non-nested changes per shell | 157 | 0 |
| Shells with one count over the whole window | 98 of 3,761 | 2,745 of 3,761 |
| Equator counts used / changes along the plane | 3,764 / 3,763 | 52 / 51 |

Tying the equator count to blocks of shells (1.5%, 3%, 6% tolerance) changed nothing for the chain, because its counts are already shared by runs of shells.

### 3.3 Equator counts [C]

`08_equator_count_sets.py`: allowed equator counts `6m` for m built from small primes, wedge counts 630 and up.

| m made of | Counts | Biggest step | Equator arc error | Median share of meridians kept across a change |
|---|---|---|---|---|
| 2 and 3 | 34 | 18.5% | up to 9.3% | 12.5% |
| 2, 3 and 5 | 100 | 11.1% | up to 5.6% | 1.3% |
| 2, 3, 5 and 7 | 198 | 7.1% | up to 3.6% | 1.6% |

Adding the prime 5 or 7 tightens the equator but lets a chain use a ratio-5 or ratio-7 step near the pole (arcs down to 0.45 and up to 2.2 there at shell 3,000), and a change keeps almost no meridians.

## 4. How this fits the ladder

The cylinder's ladder (243/81/27/9/3/1) is built from nested wedge counts so every parent wedge line is a child line, with sectors joining the block that holds their centre. On shells the same idea needs counts that divide: a chain is exactly that, with latitude and shell directions both handled by dividing counts. A block of 3 or 9 shells can use one count where the equator count is shared (most of the galaxy beyond shell 100), and its band edges are the shared latitude steps; blocks near the poles sit on the chain's coarse counts. The multiples-of-6 rule gives blocks only the six pole meridians and would need the cylinder's master-wedge trick on top (and then loses some of its shape advantage).

The pole cap and core are unchanged: the core ball is radius 2 pc (diameter 4 pc), each shell has two caps with 6 rim wedges (radius 3.82 pc), as approved.

## 5. Questions for Boss

1. **Inside a shell: fine counts or nested counts?** Multiples of 6 (arcs ±0.1%, hundreds of broken meridians per shell) or a nested chain (arcs 0.73 to 1.27, no broken meridians, 73% of shells with one count), or the middle choice (fine in the disk window, nested above it; not measured).
2. **Along the plane: how fine may the equator count be?** One count per shell (3,763 changes), `6 * 2^a * 3^b` (34 counts, equator arcs ±9%, 12.5% of meridians kept across a change), or with prime 5 (100 counts, ±5.6%, 1% kept).
3. **How much arc error is acceptable at high latitude?** The chain's 0.71 to 1.41 swing falls on cells above roughly 25 degrees, which are 11% of cells; a tighter bound needs non-nested steps.
4. **Do blocks have to cut on cell boundaries?** If yes, the chain (or the master trick) is required; if blocks may take the cells whose centre they hold, as the cylinder ladder does, the plan's rule can stay and the ladder item picks a block count per level.

## 6. What I could not measure

- The ladder itself on shells (block keys, stage cache, block statistics): not designed here.
- The middle choice in question 1, and chains that use ratios other than the smallest primes first.
- The cell count of a stair-stepped equator in the dense core shells (1 to 3 are hand-tabled in the plan).
- Any time cost: this only counts structure.

## 7. Method and files

`06_latitude_chains.py` (per-shell chains, the table in 3.1), `07_wholegalaxy_wedge_structures.py` (about 3.5 minutes; section 3.2), `08_equator_count_sets.py` (3.3), `sphgrid2.py` (cap, band and wedge layout, shared with the plan). Output in `results-ladder-latitude.txt`, `results-ladder-galaxy.txt` and `results-ladder-equator-sets.txt`. The chain picks, for each band, the allowed count nearest the wanted `2 pi R sin(colatitude) / 4 pc` in log terms, and divides the equator count by its smallest primes first so the ratio-3 step sits nearest the pole.
