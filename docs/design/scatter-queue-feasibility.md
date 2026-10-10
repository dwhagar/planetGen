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

- **Speed.** About 78 times faster at quarter-scale (47 s to 0.6 s) and about 35 times at default scale (about 12.7 minutes to 20 s). Measured later (see Layer stacks below): the majorant over all 3.1M rings is 3.75 s of the run, the rest is the candidates (about 37 µs each, 41% of them rejected), so the time is spent on candidates and objects, not on rings.
- **Counts.** 12 seeds: mass mean 4,189.6 (sd 53.4), luminosity 4,996.2 (sd 40.1), +0.7% and +0.6% against expected, within about 1.5σ.
- **Spatial distribution.** At a 1,000 Lsun floor (about 420k stars per replicate, a 9 × 8 × 8 histogram) the chi-square per degree of freedom was 0.93 to 1.02 between the sampler and today's code, 1.07 pooled (p about 0.12). Ring-band and layer-band shares agree within 0.05 points. The two cannot be told apart.

## Follow-up (Boss, 2026-10-10): several objects per sector, and stacks of layers

### Several objects per sector: tiered capacity

Boss's variant: keep several objects per sector where that is what we have today, still place by popping from a queue, and give dense sectors a higher capacity by density tier. Evaluated on the per-sector expected counts at default scale (2.39e10 sectors, arrays in `scripts/lamarr_*.npz`, analysis in `scripts/tier_analysis.py`, output in `scripts/tiers.log`).

**Capacity is a function of the sector's expected count λ.** Take the smallest capacity from 1, 2, 4, 8, 16 such that P(Poisson(λ) > capacity) is below ε. At ε = 1e-4 the tiers fill up like this:

| Pass | Sectors at capacity 1 | 2 | 4 | 8+ | Objects dropped by the cap | Queue slots per sector |
|---|---|---|---|---|---|---|
| Mass (8 Msun) and luminosity (5,000 Lsun) | 100% | | | | 0.02% (same as cap 1) | 1.000 |
| Luminosity, 1,000 Lsun | 98.6% | 1.38% | 0.019% | | 0.13% (was 0.84% at cap 1) | 1.014 |
| Luminosity, 300 Lsun | 82.0% | 14.9% | 2.55% | 0.57% | 0.07% (was 11.8%) | 1.265 |
| Luminosity, 100 Lsun | 70.1% | 23.5% | 5.06% | 1.33% | 0.05% (was 18.9%) | 1.49 |
| Phenomena (14 Msun) | 96.3% | 3.37% | 0.29% | | 0.12% (was 1.43%) | 1.043 |

So tiers fix the low-floor and phenomena loss, and cost nothing at the shipped star floors, where every sector sits in tier 1. A tier is also a clean density rating: a sector's count never exceeds its capacity, so its fill (count over capacity) stays between 0 and 1.

**The queue itself changes the statistics; independent draws with a tier cap do not.** If a sector is listed `c` times and each pop takes one slot, a sector with expected count λ ends up with Binomial(c, λ/c) objects, not Poisson(λ). That has the same mean but fewer sectors with two or more objects: with c = 2 the chance of two is λ²/4, half of Poisson's λ²/2. Measured on the default phenomena load, objects beyond the first in a sector are 1.43% under Poisson, 0.50% for a queue with ε = 1e-3 tiers, and 0.81% with ε = 1e-4 (at the 100 Lsun floor 18.9% against 17.3%). Cap 1 gives none. Two ways to avoid it:

1. **Preferred: draw as today (independent, object-first), count objects per sector in a small dictionary, and drop an object when its sector is at capacity.** The count is Poisson up to the cap, so dispersion is exactly today's, and the cap is a pure function of λ. At ε = 1e-4 the drop is 0.05% to 0.13% of objects. Cost: one dictionary lookup per accepted object.
2. A true queue needs capacity much larger than λ (copy weight λ/c with c ≥ 8λ) before it looks Poisson, which means long lists. Not worth it.

The tier is computed from λ at the sector centre, so it is the same for every object in the sector. Not measured: the sampler running with the cap (only the analytic effect above; the queue's Binomial shape is a standard result and was not simulated).

### Stacks of layers (the sparse method)

Boss's sparse method for runs of blank layers, adapted: one majorant, one Poisson draw and one ring table for a stack of G consecutive layers, layer chosen uniformly among the stack's layers that reach the sampled ring, accept with true density over the stack majorant (`scripts/stack.py`, `bench_stack.py`).

| Scale | Stack size G | Time | Objects | Candidates |
|---|---|---|---|---|
| Quarter | 1 | 1.11 s | 12,484 | 21,456 |
| Quarter | 4 | 0.93 s | 12,285 | 23,049 |
| Quarter | 16 | 1.07 s | 12,519 | 29,906 |
| Quarter | 64 | 2.12 s | 12,555 | 67,723 |
| Quarter | 509 (all) | 4.97 s | 12,442 | 150,546 |
| Default | 1 | 49.5 s | 795,968 | 1,350,353 |
| Default | 4 | 46.0 s | 795,528 | 1,374,921 |
| Default | 16 | 49.4 s | 795,313 | 1,479,856 |

(Luminosity pass at 5,000 Lsun with an 8 Msun mass split, so the object counts differ from the 14 Msun tables above. Single process.)

**Result: stacks do not pay.** A bigger stack has a looser majorant, so more candidates are rejected, and that cancels the saving of fewer majorants. Best case is G = 4 at 7% faster; G = 16 is break-even; whole-galaxy stacks are 5 times slower at quarter-scale. Reason: an empty layer costs only its majorant, about 1.2 µs per ring (3.75 s for all 3.1M default rings), so there is nothing left to skip. What costs time is the candidates, about 37 µs each (11.5 µs to place the point and check its address, 5.9 µs for the density, the rest loop overhead), with 41% rejected, and the per-object work after that (about 62 µs per object here).

So the stack idea's job, avoiding a visit to every empty layer, is done by the object-first draw itself: each layer costs one majorant and one Poisson draw whether it is empty or not. PERF.57's grouping is not needed.

**Where a further 10% to 20% is, if wanted.** Rejected candidates are 41% of the 1.35M and cost the same 37 µs as accepted ones. Doing the sector-address check only after the density test saves about 8 µs per rejected candidate (about 4.5 s of 49 s). A tighter majorant (per azimuth bin, not per ring) would raise the 59% acceptance toward 85%. Neither is needed for the 35 times.

**Recommended wording for PERF.58.** Per-layer sampler (a stack of one); stack size is a named tuning value defaulting to 1, so a stack of G can be tried later without a code change (G = 4 measured 7% faster). Several objects per sector: independent draws with a per-sector count dictionary and a tier cap from λ at the sector centre (smallest of 1, 2, 4, 8, 16 with P(N > cap) < 1e-4), objects over the cap dropped; the cap applies in the star passes and phenomena alike and replaces the one-per-sector rule. Optional: defer the address check until after acceptance.

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
5. **Use tiered capacity, not cap-of-one,** so several objects per sector stay possible (see Follow-up): independent draws, a per-sector count, cap from λ.
6. **Do not group layers into stacks.** Measured: no gain over one layer at a time (see Follow-up).

Suggested TODO items (for the TODO thread to file if Boss agrees): share `_ring_bins` across passes; object-first thinned sampler for the star passes (supersedes the PERF.57 grouping); large-mean Poisson helper; object-first sampler for phenomena.

## Sources

Scripts and logs: `research/scatter-queue/scripts/` (`profile_today.py`, `lambda_dist.py`, `bench_literal.py`, `thinned.py`, `bench_thinned.py`, `seeds_thinned.py`, `equiv.py`). Related: [generation-performance-study.md](generation-performance-study.md); TODO items PERF.57, GEN.187, GEN.195.

## Phenomena prototype result (Bugfixes lane 1, 2026-10-10, PERF.61)

The object-first sampler was built for the phenomena pass (per-class Poisson counts from the stellar density bound times a per-kind regional factor bound, acceptance on the true density times the regional factor, tiered sector capacity across classes) and measured against today's `phenomenon_scatter.scatter_layer` on the default-scale galaxy (2,041 layers, outer ring up to 3,763, 14 Msun cut, single process; scripts and log in `research/scatter-queue/scripts/bench_phenomena_object_first.*`, prototype in `phenomenon_scatter_object_first_prototype.py`).

| | Today | Object first |
|---|---|---|
| Objects (41 sampled layers) | 1,385,251 | 1,380,928 (-0.3%) |
| Time (41 sampled layers) | 41.9 s | 42.1 s |
| Whole pass, single process (scaled) | 34.8 min | 35.0 min |

- **Counts and positions match** (per-layer counts within 0.5%).
- **No speed gain.** A default-scale layer holds about 33,000 to 95,000 phenomena, so the pass is dominated by per-object work, not by the ring walk: today about 17 us an object (bin weights are precomputed), object first about 30 us (a candidate costs a point draw, a density evaluation and a regional factor, and about half of them are rejected). Only sparse outer layers gain (layer 700: 0.32 s to 0.08 s), and the dense ones lose (layer 0: 1.4 s to 2.35 s).
- Deferring the sector-address check until after acceptance and using `relative_density` instead of the population split cut the new sampler's time by about 40% (3.9 s to 2.35 s for layer 0), and it still does not beat today's.
- **Decision:** the phenomena pass keeps today's ring walk; the prototype is not merged. A further gain would need the candidate work vectorised in numpy (points, density, regional factor and acceptance for thousands of candidates at once), a larger change that would need its own prototype; not done.
- The star passes keep the sampler (PERF.58), where objects are few (267,000 to 318,000) and the ring walk was the whole cost.
