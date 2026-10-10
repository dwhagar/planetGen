# Scatter "sector queue" feasibility study

Question (Boss, 2026-10-10): would a shuffled list of every available sector, popped one object at a time, make scattering faster, and what would it change?

**Short answer.** The literal list is slower, by roughly two orders of magnitude. The idea behind it (decide how many objects there are, then pick where, instead of visiting every ring) is right, and a prototype built on it ran 35 times faster per star pass at default scale with a spatially identical result. It does not need a list, and it does not need one-object-per-sector.

Method and scripts: `/mnt/project-files/research/scatter-queue/scripts/` (all numbers below come from there). Where the existing cost model applies, this note cites [generation-performance-study.md](generation-performance-study.md) rather than redoing it.

## 1. Where scatter time goes today

Quarter-scale galaxy (509 layers, outer ring 940, 3.74e8 sectors), single process, warm, no database, measured per pass:

| Pass | Wall | Objects | Time in `_ring_bins` | Other |
|---|---|---|---|---|
| Mass | 46.9 s | 4,083 | 43.6 s (93%) | 461 of 509 layers empty, 35.7 s spent in those |
| Luminosity | 47.5 s | 4,935 | 43.9 s (92%) | 389 empty layers, 22.8 s |
| Phenomena | 73.8 s | 1,375,520 | 42.7 s (58%) | place 16.0 s, kind weights 8.5 s, Poisson 0.9 s |

- The cost is the ring walk, not the objects. Each ring evaluates the density model about 32 times (about 6 µs per call, 0.22 to 0.24 ms per ring). The passes find a few thousand stars and still pay for every ring of every layer.
- All three passes recompute identical `_ring_bins` inputs. Sharing them would save about 87 of 168 s at quarter-scale.
- At default scale (2,041 layers, 3,136,126 ring visits per pass) the measured 0.244 ms per ring is about **12.7 minutes per pass** for the two star passes.
- The per-object costs from the performance study (draw 40 to 49 µs, insert 38 µs) only matter for the phenomena pass, which places millions of rows. No database writes were measured here; the study's insert cost applies unchanged to any design.
- `_sample_poisson_count` is Knuth's O(mean) algorithm. It is fine per ring, but any design that draws one count per layer needs a large-mean Poisson helper.

## 2. Cost of the list approach at default scale

Default galaxy: 2.3921e10 candidate sectors.

- **Building and shuffling.** Layer 0 at quarter-scale (2.78M sectors): a Python tuple list takes 2.07 s and 267 MB, shuffle 1.75 s. A numpy int32 index array takes 11 MB and shuffles in 0.27 s. The largest default layer (4.4e7 sectors) is 4.2 GB as tuples, 178 MB as int32. So the list can exist if it is an index array, not tuples.
- **Weights.** Step 5 needs each sector's density. At 6.07 µs per sector, layer 0 costs **16.9 s** against 0.246 s for today's ring walk (about 68 times more). Across the default galaxy that is 2.39e10 × 6 µs, roughly **40 CPU-hours per pass**, against 12.7 minutes today (190 to 270 times more density evaluations).
- **Weighted shuffle.** Efraimidis-Spirakis keys plus `argpartition` cost 0.28 s per layer once the weights exist (peak 528 MB), so the shuffle is cheap and the weights are not.
- **Sparsity makes the list unnecessary.** The mass and luminosity passes place about 2.7e5 and 3.2e5 objects among 2.4e10 sectors. Picking k distinct sectors needs a used-set or a lazy Fisher-Yates over k draws, not N keys. The list buys nothing over sampling.

## 3. What would change in the results

Expected objects per sector (λ) at default scale:

| Pass | Expected objects | Sectors with λ ≥ 0.001 | Lost under one-per-sector cap |
|---|---|---|---|
| Mass (8 Msun) | 266,390 | 0.12% | 0.02% |
| Luminosity (5,000 Lsun) | 317,800 | 0.10% | 0.02% |
| Phenomena (14 Msun cut) | 6.7155e7 | 4.57% have λ ≥ 0.01; 0.146% have λ ≥ 0.1 (6.1% of objects) | 1.43% |

- **One object per sector costs nothing for stars.** The largest per-sector λ is 0.004 (mass) and 0.003 (luminosity) at the defaults, so a second object in a sector is a rare event and the cap drops 0.02% of objects.
- **It matters for phenomena (1.4% lost) and for low luminosity floors.** With the luminosity floor at 1,000 Lsun the cap loses 0.84% of 2.6e7 objects; at 100 Lsun it loses 18.9% (largest λ 3.3). Capping at one makes those regions come out sparser than the physical density, and the "move it to another sector" variant over-counts density next to filled regions. If the cap is adopted it should be applied only to the star passes at the shipped floors, or the sampler should draw the number per sector from Poisson(λ) rather than 0 or 1.
- **Dropped objects.** Today an object that lands in a filled or mass-marked sector is drawn and dropped. That equals thinning when the expected count is computed over the available sectors only, so the dropped-object behaviour matches Boss's step 4 (do not put those sectors in the list) as long as the expected count is not recomputed over the removed sectors and renormalised.
- **Seeds.** Any new sampler is a different random sequence, so existing seeds give a different galaxy. A reseed (`planetgen plan`) is needed, as with the other scatter changes.

## 4. How it fits PERF.57 and GEN.187

- **PERF.57** (grouping layers 10/20/40 after 5 empty layers, halving after a placement) is a majorant over a group of layers: it exists to avoid paying for empty layers. The object-first sampler is the same trick taken to its end, with one majorant per layer drawn once and no group size to tune. If PERF.57 lands first, the sampler replaces its group logic; if not, PERF.57 can be dropped in favour of it. Bugfixes lane 1 should know before it finishes PERF.57.
- **GEN.187** (bright-star back scatter by mass rings, 1, 2, 5, 8 Msun around filled space) works on a small explicit candidate set (face-adjacent rings of filled sectors). A per-candidate draw is natural there, and Boss's queue idea fits it well. This was checked against the TODO text only, since GEN.187 is not built.
- Boss's sparse method for runs of blank layers is the same idea at layer scale: do not visit what cannot produce anything.

## 5. Prototype: object-first thinned sampler

Per layer: draw N from Poisson(Σ e·slots·M) where M is a certified density majorant per ring; choose the ring from a cumulative table; choose slot and point as today; accept with probability true_density / M; skip filled sectors. The majorant is exact enough to be cheap: density falls with radius and |z|, the arm term is bounded by 1 + A, and each population's share peaks at cos = ±1. A `check=True` assertion that true density never exceeds M was never violated.

| Scale | Pass | Time | Objects | Expected |
|---|---|---|---|---|
| Quarter | Mass | 0.61 s | 4,298 | 4,161 |
| Quarter | Luminosity | 0.69 s | 5,116 | 4,964 |
| Default | Mass | 19.6 s | 267,092 | 266,390 |
| Default | Luminosity | 21.0 s | 318,218 | 317,800 |

- **Speed.** About 78 times faster at quarter-scale (47 s to 0.6 s) and about 35 times at default scale (about 12.7 minutes to 20 s). The default time is mostly the numpy majorant over 3.1M rings, which could be cached across passes and runs.
- **Counts.** 12 seeds: mass mean 4,189.6 (sd 53.4), luminosity 4,996.2 (sd 40.1), +0.7% and +0.6% against expected, within about 1.5σ.
- **Spatial distribution.** At a 1,000 Lsun floor (about 420k stars per replicate, a 9 × 8 × 8 histogram) the chi-square per degree of freedom was 0.93 to 1.02 between the sampler and today's code, 1.07 pooled (p about 0.12). Ring-band and layer-band shares agree within 0.05 points. The two cannot be told apart.

## What could not be measured

- No database writes: the study builds the galaxy shape in memory, so insert and commit costs are excluded (see the performance study for those).
- Single process on a 4-core container; the production worker pool and PERF.42 pre-warm were not run.
- The phenomena pass was not prototyped. It needs per-kind majorants and still pays placement and row costs (about 16 s per 1.4M objects at quarter-scale).
- The prototype uses numpy's Poisson for the layer count and evaluates density at the exact point rather than at bin centres.
- GEN.187 and PERF.57 are not on main, so the interaction is analysis, not measurement.

## Recommendation

1. **Do not build the shuffled sector list.** It costs 190 to 270 times more density evaluations and gains nothing, because the objects are sparse.
2. **Share `_ring_bins` across the three passes.** Low-risk, no change to results, about a 2 times saving at quarter-scale.
3. **Build the object-first sampler for the mass and luminosity passes**, replacing the per-ring walk and PERF.57's grouping. About 35 times faster at default scale, distribution unchanged, needs a reseed and a large-mean Poisson helper.
4. **Phenomena pass next, separately.** Same idea with per-kind majorants. The gain is smaller because row costs dominate.
5. **Keep the cap-of-one only where λ is small** (the shipped floors); do not apply it at low luminosity floors or to phenomena.

Suggested TODO items (for the TODO thread to file if Boss agrees): share `_ring_bins` across passes; object-first thinned sampler for the star passes (supersedes the PERF.57 grouping); large-mean Poisson helper; object-first sampler for phenomena.

## Sources

Scripts and logs: `research/scatter-queue/scripts/` (`profile_today.py`, `lambda_dist.py`, `bench_literal.py`, `thinned.py`, `bench_thinned.py`, `seeds_thinned.py`, `equiv.py`). Related: [generation-performance-study.md](generation-performance-study.md); TODO items PERF.57, GEN.187, GEN.195.
