# Nebula density against the real Milky Way, and what to generate versus what to draw

Boss (2026-10-10 22:01Z): "review and have a thread look at the cloud density vs the real world and how we can manage that cloud density for making things visible."

Follows [nebula-map-visibility.md](nebula-map-visibility.md) (the three-scale map design). Informs GEN.152 (lower the cloud field to about 1%), DB.24 (the `nebula_field` table), MAP.172 to MAP.178 (aggregates, tile layer, regions, cover, field clouds), GEN.202 (scattered nebulae), MAP.182 (hosted H II for pre-placed stars). Builds on [nebula-and-asteroid-field-classes.md](nebula-and-asteroid-field-classes.md) (the GEN.47 field and its first rate check).

Status: research, 2026-10-10; nothing here is built. Evidence tags: [S] seen in the repo or in a cited source, [C] computed or measured in this research, [R] recalled and not re-checked (listed at the end).

## Summary

- **The field is too dense in counts, about right in packaging, and the whole fix is a handful of constants.** Against the Milky Way, class M (giant molecular clouds) is about 10 times too many (102,000 against 8,107 to 9,710 catalogued), class N (dark clouds) is in range by count (about half the Lynds catalogue's entries within 1 kpc) but its clouds are large, and classes P and Q (Bok globules, cores) are in the right range, and the sphere volume the field fills is 2% / 8% / 19% / 24% where the gas factor is 0.2 / 0.9 / 2.3 / 3.0, against 0.5% to 1% observed for the thin disc and 1% to 2% in arms [C, S]. (This corrects the first rate check, which said 3 / 13 / 35 / 46%: it gave the four classes equal weight, but the field draws them 1:3:2:1.)
- **Recommended generation rates, all in the GEN.152 "about 1%" default Boss approved:** keep 9.5% of class M (9,700 galaxy-wide at the default scale), keep 50% of class N (153,000, a volume decision, see section 6), keep P and Q whole (192,000 and 99,000). The field then fills 0.2% / 1.1% / 2.8% / 3.6% of volume (bounding spheres) and holds 454,000 clouds instead of 700,000. Applying it as an acceptance test on the existing draw, so surviving clouds keep their cells and object IDs, costs one hash per cloud.
- **Storing the field costs far less if only classes M and N are in the table.** P and Q make 291,000 of the 454,000 clouds, are 0.2 to 1.6 pc across (invisible beyond about 250 pc), and are cheap to derive from the seed when a sector or a near view needs them. `nebula_field` then holds about 163,000 rows (about 11 MB) instead of 700,000 (about 50 MB).
- **Lowering density does not shrink the map's tiles; the budget does.** In a dense 4 kpc sample, regions at 16 px number 3,142 at 8 kpc in the current field and 2,917 after the cut, against a draw budget of 600 (MAP.176). Tile size is set by that budget and by the size law. The cut matters for the database, for how many sectors sit inside a cloud, and for how much sprite overlap the near views carry.
- **Draw side: three small rules.** (1) A far sprite is drawn at the cloud's equivalent radius (0.49 of its bounding radius, because the shape fills only 12% of its bounding sphere [C]), with the existing 2 px floor; the bounding radius stays for picking. (2) Region budget 600, ranked by volume. (3) The Nebula cover ramp is re-binned to the realistic density (0, under 0.5%, 0.5 to 1.5%, 1.5 to 3%, 3% and up), because the first draft's 5 / 15 / 30% bins would leave almost every block in the lowest two.
- **Two other kinds are also off.** Supernova remnants in the scatter (1e-8 per pc^3) number about 20,000 against about 300 catalogued and about 400 to 3,000 expected from the supernova rate; 1.5e-9 gives about 3,000. Hosted H II regions follow the O-star count (every O star gets one): at 270,000 O stars they would number 270,000 against about 8,000 catalogued; a chance of 0.25 per O star gives the catalogue's ratio. Planetary nebulae (28,000) are inside the published range of 4,000 to 46,000 and stay.
- **This changes Boss's generated galaxy, flagged:** a new galaxy takes the new rates from the plan. An existing one keeps its stored clouds until a one-off cleanup removes the 35% of stored field clouds the new rule would not have made; or he re-plans, as he does for timing runs. Details in section 8.
- Eleven build items, with phases, are in section 9. Defaults in section 10.

## 1. What the generator makes now

[S] `nebula_field.py`, `tuning.py`; [C] this research.

The galaxy is cut into 50 pc cells. Each cell draws a Poisson number of dark-family clouds with mean `PHENOMENON_DENSITY_PC3["molecular-cloud"]` (5e-6 per pc^3) times the cell volume times the gas factor (the young-star tracer with the gas's own arm contrast, to the power 1.4, capped at 3). Each cloud is a `Nebula(nebula_type="dark")`, whose class is drawn with the `NEBULA_CLASSES` frequencies 1 : 3 : 2 : 1.

| Class | Name | Share of draws | Radius, ly (p10 / median / p90) | Median radius, pc | Sphere volume, pc^3 (mean) | Galaxy-wide now [C] |
|---|---|---|---|---|---|---|
| M | Giant molecular cloud | 14.6% | 35 / 81 / 137 | 25 | 120,000 | about 102,000 |
| N | Dark cloud | 43.8% | 5.5 / 26 / 45 | 8 | 3,900 | about 307,000 |
| P | Bok globule | 27.5% | 0.6 / 1.6 / 2.7 | 0.5 | 1 | about 192,000 |
| Q | Star-forming core | 14.1% | 0.2 / 0.5 / 0.9 | 0.2 | 0.1 | about 99,000 |
| | | | | | Total | about 700,000 |

Two properties matter for everything below:

1. **Only M and N carry volume.** Class M is 14% of the clouds and 91% of the filled volume at equal gas factor; M and N together are above 99.9%.
2. **The shape fills 12% of the bounding sphere.** `radius_ly` is the distance to the shape's farthest point. Sampling 500 generated clouds, the metaball shape fills 12.2% of that sphere on average (class M 12.2%, N 12.7%, P 12.9%, Q 12.3%; p10 8%, p90 18%) [C, `shape_fill_fraction.json`]. So a cloud's true volume is about an eighth of its bounding sphere, and the sphere is what a sprite draws and what the 2 px size law tests.

Filling by environment (share of volume inside bounding spheres, `1 - exp(-sum of n V)`) [C, the gas factors are the first rate check's 0.2 / 0.9 / 2.3 / 3.0 at the solar circle]:

| Gas factor | 0.2 (between arms) | 0.9 (midway) | 2.3 (arm crest) | 3.0 (cap: inner disc) |
|---|---|---|---|---|
| Now | 1.8% | 8.0% | 19.2% | 24.2% |
| Sampled box 6 kpc from the core (mean over the box) | 11.0% | | | |

The old note's 3 / 13 / 35 / 46% used a mean cloud volume of 3.0e4 pc^3, which weights the four classes equally. The field draws 1 : 3 : 2 : 1, so the mean is 1.7e4 and the figures above are right.

## 2. What the real Milky Way has

All counts are galaxy-wide unless stated; links are the sources I could confirm in this research [S] unless tagged.

| Object | Real count and size | Source |
|---|---|---|
| Giant molecular clouds | 8,107 clouds holding 98% of the CO emission, total H2 mass about 1.2e9 Msun, most probable radius about 30 pc; the author's later release lists 9,710 | Miville-Deschenes, Murray and Lee 2017, [arXiv 1610.05918](https://arxiv.org/pdf/1610.05918), [catalogue page](https://dc.g-vo.org/rr/q/lp/custom/CDS.VizieR/J/ApJ/834/57) |
| Massive GMCs | 1,064 | Rice et al. 2016, quoted in the GEN.47 note [R] |
| Molecular gas fraction in bound clouds | 40% of CO mass | Heyer and Dame 2015, via [arXiv 1812.02180](https://arxiv.org/pdf/1812.02180) |
| Cold medium volume filling, thin disc | under 5% (a simulation cross-check, not a measurement) | [arXiv astro-ph/0407034](https://arxiv.org/pdf/astro-ph/0407034) |
| Molecular volume filling | about 0.5% inner Galaxy, about 1% thin disc; arms 1% to 2% | the GEN.47 note [S/C] and `GMC_ARM_FILLING_FACTOR` = 0.015 in `tuning.py` [S] |
| Dark nebulae (Lynds) | 1,802 entries, mostly within about 1 to 2 kpc [R for the range] | [Lynds 1962](https://articles.adsabs.harvard.edu/pdf/1962ApJS....7....1L), [catalogue](https://en.wikipedia.org/wiki/Lynds%27_Catalogue_of_Dark_Nebulae) |
| Bok globules | 248 optically selected northern globules catalogued, 169 southern; a Galactic total of 3.2e5 small globules (a secondary citation of Clemens et al. 1991, so [R]) | [Bok globules in the LMC, arXiv astro-ph/9902165](https://arxiv.org/pdf/astro-ph/9902165) |
| Cold clumps (cores) | 13,188 Planck sources, a survey-limited floor | Planck 2015 XXVIII, [arXiv 1502.01599](https://arxiv.org/pdf/1502.01599) |
| H II regions | over 8,000 in the WISE catalogue, about 1,500 confirmed by recombination lines or H-alpha, about 2,500 more with radio continuum; a later census counts about 7,000 with 2,000 confirmed | Anderson et al. 2014, [arXiv 1510.07347](https://arxiv.org/pdf/1510.07347), [catalogue](https://dc.g-vo.org/rr/q/lp/custom/CDS.VizieR/J/ApJS/212/1) |
| O stars | 14,000 to 18,000 O systems in a basic model, 30,000 to 50,000 individual O stars allowing for binaries and extinction; 24,706 OB stars within 1 kpc | Maiz Apellaniz et al. 2013 and Quintana, Wright and Garcia 2025, [arXiv 2503.08286](https://arxiv.org/pdf/2503.08286) |
| Core-collapse supernova rate | 0.4 to 0.5 per century (Quintana et al.); 3 per century is the older rate used in the CTA planning, [arXiv 2310.02828](https://arxiv.org/pdf/2310.02828) | |
| Supernova remnants | 294 in the revised Green catalogue, 310 on its current page; the count expected from the supernova rate over a remnant life of about 1e5 years is 400 (0.4 per century) to 3,000 (3 per century) [C] | [Green 2019, arXiv 1907.02638](https://arxiv.org/pdf/1907.02638), [catalogue](https://www.mrao.cam.ac.uk/surveys/snrs-2022/) |
| Planetary nebulae | about 3,000 known; estimates of the total run 4,000 to 46,000: 7,200 plus or minus 1,800 (Peimbert), 17,000 plus or minus 3,000 (Bobylev and Bajkova), 46,000 plus or minus 13,000 (Moe and De Marco), about 25,000 typical | [arXiv 1704.01718](https://arxiv.org/abs/1704.01718), [arXiv 0910.0465](https://arxiv.org/pdf/0910.0465) |

Counts of clouds are sharper than filling fractions (a filling needs a layer volume and a shape for filamentary gas), so the counts are the anchor and filling is the cross-check.

## 3. Ours against the real one

| Kind | Ours, default scale | Real | Factor | Verdict |
|---|---|---|---|---|
| Class M | about 102,000 | 8,107 to 9,710 | about 10 to 12 | too many |
| Class N | about 307,000 | no clean total; Lynds 1,802 over the nearby sky. Our N density at gas 0.9 is 2e-6 per pc^3, which puts about 900 N clouds in a 1 kpc radius, 150 pc thick slab round the Sun against Lynds' 1,802 entries, mostly smaller clouds [C, R] | about 0.5 | count in range or a little low; the clouds are larger (median 8 pc against mostly 0.5 to 5 pc) and add 0.75% of the 1.5% filling at midway |
| Class P | about 192,000 | 3.2e5 (secondary citation) | 0.6 | right range |
| Class Q | about 99,000 | 13,188 Planck clumps is a floor; dense cores are far more numerous | | right range, unconstrained |
| Sphere volume filled, midway / crest | 8.0% / 19.2% | 0.5% to 1% disc, 1% to 2% arms | about 8 to 10 | too full |
| H II regions (hosted: one per O star, half of B0 to B2) | about 270,000 O stars at the default scale [R scaling of 4,217 at quarter scale] gives about 270,000 plus B0 to B2 hosts | about 8,000 catalogued | about 30 | too many, because it follows the O-star count |
| Supernova remnants (scatter) | about 20,000 | 310 catalogued, 400 to 3,000 expected | 7 to 50 | too many |
| Planetary nebulae (scatter) | about 28,000 | 4,000 to 46,000 | 0.6 to 7 | inside the range |

The shape factor is why counts and volume do not agree in the real direction. At the catalogue count of M, a mean true cloud volume of 1.2e5 pc^3 x 0.12 = 1.45e4 pc^3 would fill only 0.13% at gas 0.9 [C], below the observed 0.5% to 1%. Real clouds are therefore either larger than ours at equal count or the observed filling is for a thinner layer. I cannot resolve that from the sources I have, so the recommendation anchors counts (well measured) and the sphere filling the project already targets (GEN.152's "about 1%"), and lists a radius increase as an option, not a default (section 4, option F).

## 4. The options, measured

A 4 kpc box on the plane centred 6 kpc from the core holds 45,061 clouds (dense: it spans the solar circle and the inner arm; the mean gas factor is about 1.4). Poisson thinning by a per-class keep probability is exactly a Poisson draw at the lower rate, so each option thins that one sample [C, `density_options.py`].

| Option | Keep M / N / P / Q | Galaxy-wide clouds | Filling by bounding sphere at gas 0.2 / 0.9 / 2.3 / 3.0 | Sampled box filling |
|---|---|---|---|---|
| A Now | 1 / 1 / 1 / 1 | 700,000 | 1.8 / 8.0 / 19.2 / 24.2% | 11.0% |
| B Uniform x0.1 | 0.1 each | 70,000 | 0.2 / 0.8 / 2.1 / 2.7% | 1.2% |
| C M to the catalogue | 0.095 / 1 / 1 / 1 | 608,000 | 0.3 / 1.5 / 3.7 / 4.8% | 2.2% |
| **D M and N cut (recommended)** | **0.095 / 0.5 / 1 / 1** | **454,000** | **0.2 / 1.1 / 2.8 / 3.6%** | **1.7%** |
| F D plus class M radius x2 | 0.095 / 0.5 / 1 / 1 | 454,000 | 1.4 / 6.0 / 14.5 / 18.5% | 8.6% |

What a map draws, for the sample box (individual clouds with a bounding radius of at least 2 px; regions are groups of clouds in one 16 px cell that span at least 2 px; "ink" is summed sprite area over a 900 x 1400 px screen, so above 1 means overlaps; "eq" scales each radius by 0.49) [C]:

| Camera radius | A individual / regions / ink | D individual / regions / ink / eq | F individual / ink / eq |
|---|---|---|---|
| 15 kpc | 2,425 / 574 / 0.04 | 226 / 548 / 0.00 / 0.00 | 474 / 0.02 / 0.01 |
| 8 kpc | 4,683 / 3,238 / 0.18 | 457 / 2,917 / 0.02 / 0.00 | 586 / 0.07 / 0.02 |
| 4 kpc | 15,049 / 10,981 / 0.94 | 5,217 / 7,780 / 0.17 / 0.04 | 5,246 / 0.39 / 0.09 |
| 2 kpc | 20,700 / 6,122 / 3.91 | 7,884 / 2,929 / 0.74 / 0.18 | 7,884 / 1.60 / 0.38 |
| 1 kpc | 14,404 / 685 / 9.88 | 5,656 / 309 / 1.84 / 0.44 | 5,665 / 4.10 / 0.98 |
| 500 pc | 4,416 / 25 / 11.70 | 1,781 / 12 / 2.30 / 0.55 | 1,787 / 5.47 / 1.31 |
| 250 pc | 1,405 / 1 / 13.01 | 704 / 1 / 2.25 / 0.54 | 709 / 6.44 / 1.54 |

Reading it:

- **A (now) smothers the near and middle views.** At 1 kpc and below a screen carries 10 to 13 layers of cloud sprite; between 4 and 1 kpc the 200 per tile cap leaves 800 to 1,600 of the 14,000 to 21,000 clouds that are above 2 px.
- **B is the simplest but cuts P and Q by 90%,** though both are in the real range, and takes N below its catalogue count.
- **C leaves the total heavy in arms:** 1.5% at midway but 3.7% at an arm crest, 2.5 times the 1.5% arm target. Class M alone is 1.8% at a crest, so the arm target cannot be met by cutting N: M at the catalogue count is already above it.
- **D matches both anchors:** M at the catalogue count, the total near 1% at midway, 2.8% at a crest (spheres; the shape fills an eighth of that).
- **F matches the true-volume reading but puts the spheres back at 6% to 18%,** so it undoes the visibility gain; keep it as an option only if Boss wants the true volume at 0.5% to 1%.
- **At 15 kpc and above only M reaches 2 px** (226 individual clouds in the 4 kpc box, the larger M). A realistic galaxy view therefore needs the cover shading and the regions to carry the gas; individual sprites cannot, in any option.
- **Region counts stay above the 600 budget from 8 kpc to 1 kpc in D** (2,917 to 309 at 16 px), so the budget and the volume rank of MAP.176 are required whatever density is chosen.

## 5. Generate less, or draw less?

Both, on different parts.

**Lower what is generated** where the count is wrong and the object matters in play:

- **Class M and N.** A cloud is not decoration: a sector inside one is in a region of 10 to 100 magnitudes of extinction and 10^2 to 10^6 H2 per cm^3 (the class table). At the current field 13% of arm sectors and 4% of between-arm sectors are inside a bounding sphere [S, GEN.47 note, using the old mean]; the cut lowers the filling to 0.14 of its old value, so those become about 2% and 0.6% [C scaling, not measured per sector]. Cutting M and N at the source is the only lever that changes it.
- **Supernova remnants and hosted H II regions.** Both are objects with stars, ages and pages; 20,000 remnants and 270,000 H II regions swamp lists and the Sector Map. They are rate constants (section 6).

**Keep the data realistic and control only the drawing** where the object is small, numerous and not needed as a row:

- **Classes P and Q** (291,000 clouds, 0.2 to 1.6 pc across): they stay in the seeded field exactly as now, are not put in the map's index, and are drawn only where they can be seen (near views and the sector view). The index and pyramid then hold M and N only.
- **Everything that the 2 px law hides at a given zoom:** that is the size law, unchanged, and the region grouping that follows it.

**What fraction of the picture should look nebular** (D, from the sample box above; Boss asked for "what fraction at each scale"):

| Scale (camera radius) | What carries the gas | Target on screen |
|---|---|---|
| Galaxy view, 15 kpc and out | Nebula cover shading; at most 600 regions | no individual sprites; cover top bin (3% and up) only on arm crests and the inner disc, so under 10% of the lit disc |
| Middle, 8 kpc to 1 kpc | regions and the biggest clouds | sprite ink at the equivalent radius 0.00 to 0.44 (a quarter to a half of a layer deep in arms), 300 to 600 regions |
| Near, 500 pc and in | individual clouds, shapes from `NEBULA_MESH_MIN_PX` | physical: 1% to 3% of volume, ink about 0.5 layers; P and Q appear only below about 250 pc |

## 6. Recommended rates

| Constant | Now | Recommended | Result at the default scale |
|---|---|---|---|
| Keep fraction, class M | 1 | 0.095 | 9,700 (catalogue 8,107 to 9,710) |
| Keep fraction, class N | 1 | 0.5 | 153,000 |
| Keep fraction, classes P, Q | 1 | 1 | 192,000 / 99,000 |
| `PHENOMENON_DENSITY_PC3["supernova-remnant"]` | 1e-8 | 1.5e-9 | about 3,000 (Green 310, rate 400 to 3,000) |
| `PHENOMENON_DENSITY_PC3["planetary-nebula"]` | 1.4e-8 | keep | about 28,000 (range 4,000 to 46,000) |
| `NEBULA_HOST_RULES` chance, O stars | 1.0 | 0.25 | about 1 H II per 4 O stars; with 30,000 O stars that is 7,500 against 8,000 catalogued |
| `NEBULA_HOST_RULES` chance, B0 to B2 | 0.5 | 0.1 | |
| `GMC_ARM_FILLING_FACTOR` | 0.015, read by nothing | read by a test (section 7) | arm crest 2.8% by bounding sphere |

How the keep fraction is applied: `cell_clouds` keeps drawing its full list exactly as now (so each cloud's cell, index and cloud-object-ID rank are untouched), then each cloud is accepted when a hash of `(galaxy seed, "nebula-keep", cell, index)` read as a fraction is below its class's keep value. A removed cloud still counts in the ID rank, so every survivor keeps its ID. Poisson thinning by a deterministic hash gives the same distribution as a lower rate. The cost is one hash per cloud (the `Nebula()` draw is 0.013 ms [C], so dropping the draws instead would save about 35% of a small number); a re-plan can instead cut the density constant directly and give every cloud a new ID, which is simpler code and fine for a galaxy with nothing worth keeping.

Why 0.5 for N: the Lynds comparison says our N count is in range or a little low, so this cut rests on volume and size, not count. N contributes 0.75% of the 1.5% filling at midway and its median cloud (8 pc) is larger than most Lynds clouds. Keeping N whole leaves the arm crest at 3.7%, halving it at 2.8%, dropping it at 1.8% (class M alone). 0.5 is the middle, with the lowest confidence of the three numbers, and is a constant to retune.

## 7. A test that keeps the rates honest

Rates were set by research and have drifted twice (the first field at 5e-6, the first rate check's equal-weight error). A test should state the bands and fail when a constant moves out of them:

| Check | Band |
|---|---|
| Expected class M clouds, default galaxy (the integral in `field_counts.py`) | 5,000 to 20,000 |
| Bounding-sphere filling of all field clouds at an arm crest (g = 2.3) | 1% to 4% |
| Expected planetary nebulae | 4,000 to 46,000 |
| Expected supernova remnants | 400 to 5,000 |
| Hosted H II per O star (from `NEBULA_HOST_RULES`) | 0.1 to 0.5 |
| P and Q totals | 5e4 to 4e5 |

It reads only the tuning constants and the gas model, so it runs in seconds. The bands are the real-world column of section 2, widened for the catalogues' incompleteness.

## 8. Costs and the effect on existing galaxies

| Measure | Now | Recommended (D, M and N indexed) |
|---|---|---|
| Field clouds, default galaxy | 700,000 | 454,000 (35% fewer) |
| Rows in `nebula_field` (DB.24) | 700,000 (about 50 MB at 70 bytes a row [C estimate]) | about 163,000 (about 11 MB) |
| Stored `nebulae` rows per generated sector | the number of clouds that reach it | 35% fewer; each stored nebula is about 2 KB with its text columns, shape fields, ball rows and seven indexes [C estimate from `models.py` column widths, not a measured table] |
| Sectors inside a cloud, arms / between arms | about 13% / 4% | about 2% / 0.6% [C scaling] |
| Map tile size | bounded by the 600-region budget: about 72 KB of JSON or 9 KB packed for 600 regions, 3,300 groups would be about 400 KB of JSON at 8 kpc without it | the same; the density cut does not change a budgeted tile |
| Cover shading signal | median block cover 11%, p90 17% (256 pc) | median about 2%, so the ramp is re-binned (0, under 0.5%, 0.5 to 1.5%, 1.5 to 3%, 3%+) |
| Aggregates, pyramid | counts up to 700,000 | M and N only: under 200,000 clouds into the same cell levels |
| Generation CPU | `Nebula()` 0.013 ms a draw; 0.07 to 0.55 ms a 50 pc cell; `clouds_reaching` for a sector 0 to 3 ms [C] | unchanged (all clouds are still drawn; one hash more) |

**Existing galaxies.** Stored dark-family nebulae in a galaxy made before the change came from the old field. With the hash rule, the clouds that pass the new test are the clouds that were always there and keep their IDs, so nothing breaks. About 35% of stored field clouds (14.6% x 0.905 + 43.8% x 0.5, with P and Q kept) would not be in the new field. Options, in order of simplicity:

1. **Re-plan** the galaxy (Boss re-plans for timing runs; the memory note says he resets manually). Nothing to migrate. Default for a test galaxy.
2. **A one-off cleanup command** (in the family of `planetgen plan --redo-scatters`, GEN.196): for each stored field-origin dark nebula whose hash fails the new test and that holds no object of its own (no embedded stars or remnants recorded against it), delete the row, clear the `inside_nebula_id` pointers (the foreign key is already `ON DELETE SET NULL`), and refresh containment in the sectors it reached. A dry run prints the counts first. This is the path if he wants to keep a galaxy.
3. **Do nothing:** stored clouds stay as an over-dense local patch; the map and the field index do not break (the index flags stored rows), but the patch stays too full where he has already generated.

I recommend 1 for a test galaxy and 2 once there is a galaxy worth keeping; the build items cover both. The hosted H II change and the supernova-remnant change apply to systems and scatter made after them; a `Redo scatters` run (GEN.196) rebuilds the scatter rows.

## 9. Build items (phase and order)

Numbers are placeholders; the TODO thread assigns IDs. GEN.152 and DB.24 already exist and are refined, not duplicated.

**Phase 1: foundations (data, rates, storage)**

1. **GEN.152, concrete numbers.** Per-class keep fractions M 0.095, N 0.5 (lowest confidence), P and Q 1 as an acceptance hash on top of the full draw (object IDs unchanged), in `tuning.py` as `NEBULA_FIELD_KEEP`; replaces "lower it to about 1%". Done when the sample box fills 1.7% and the galaxy holds about 454,000 field clouds. Owner: Foundations lane 1.
2. **DB.24, index M and N only.** `nebula_field` holds classes M and N (about 163,000 rows); P and Q stay seeded and are drawn from the field only in near views and in sectors. Update DB.24's size and its equality test (rows equal `clouds_reaching` for M and N). Owner: Foundations lane 1.
3. **Nebula-rate audit test.** The bands in section 7 as a unit test, plus a CLI line in `planetgen plan` output. Owner: Foundations lane 1.
4. **Supernova remnants at 1.5e-9; keep planetary nebulae at 1.4e-8.** Update the comments in `tuning.py` with the sources in section 2 (the 20,000 planetary nebulae figure is already there). Takes effect on the next `plan` or `Redo scatters`. Owner: Foundations lane 1.
5. **Hosted H II chance.** `NEBULA_HOST_RULES`: O stars 0.25, B0 to B2 0.1, later B and A unchanged; a test of the per-O-star ratio. Interacts with the O/B star counts (the mass cut, GEN.185/194): the rule is a ratio, so it stays right when the counts change. Owner: Foundations lane 1.
6. **Existing-galaxy cleanup.** `planetgen plan --redo-nebula-field` (or a sibling command): dry run, then delete stored field-origin dark nebulae the new rule rejects and that hold nothing of their own, clear `inside_nebula_id`, refresh containment. A test on a small galaxy built at the old rate. Owner: Foundations lane 1.
7. **Correct the GEN.47 note** (2 / 8 / 19 / 24%, the 12% shape fill) in `nebula-and-asteroid-field-classes.md`; done in this PR as a short pointer.

**Phase 2: the visible feature (draw side)**

8. **Far sprite at the equivalent radius.** A field cloud's far sprite uses 0.49 of its bounding radius, with the existing 2 px floor and the MAP.155 1 to 4 px opacity ramp; the bounding radius stays for picking and for the point where the shape mesh loads. Belongs in MAP.178. Owner: Foundations lane 2.
9. **Nebula cover ramp at the realistic density.** MAP.177's bins become 0, under 0.5%, 0.5 to 1.5%, 1.5 to 3%, 3% and up (optionally rescaled to the view's own percentiles), with the legend line unchanged. Owner: Foundations lane 2.
10. **Region budget stress test.** A browser or performance test that draws 3,000 candidate regions and keeps 600 by volume at 60 frames a second on the reference phone (MAP.176). Owner: Foundations lane 2.
11. **Near-view P and Q.** Draw classes P and Q from the seeded field only below about 250 pc and in the sector view, per-tile cached, with a per-tile cell budget (a 500 pc cube is 1,000 cells, 0.1 to 0.5 s uncached [C]). Default: not built for the first release; the sector view shows stored P and Q already. Owner: Foundations lane 2.

## 10. Defaults and what is flagged

- **Flagged: changes Boss's generated galaxy.** The cut removes about 35% of field clouds (91% of class M, 50% of class N) from every new galaxy; a re-plan or the cleanup command (item 6) applies it to an existing one. GEN.152 already carries Boss's approved default ("lower it to about 1%"), and D is that default made specific; nothing here goes beyond it except the class split and the two scatter and hosting constants.
- **Keep fractions M 0.095 / N 0.5 / P 1 / Q 1.** Constants; tune with the audit test.
- **Counts are the anchor, not filling.** The true-volume reading (shape fill 12%) says the cut leaves the shape volume below the observed 0.5% to 1% by a factor of 4 to 8; option F (class M radius x2) would match both but returns the bounding spheres to 6% to 18%. Not recommended.
- **Index M and N only; P and Q seeded and implicit.**
- **Regions: 16 px cells, 600 budget** (unchanged from nebula-map-visibility.md).
- **Supernova remnants 1.5e-9, planetary nebulae unchanged, hosted H II chance 0.25 (O) and 0.1 (B0 to B2).**
- **No Boss question needed;** every item has a default.

## Limits of this study

- Counts come from the in-memory field model and one 4 kpc sample box 6 kpc from the core (a dense case, mean gas factor about 1.4); outer arms and the inner disc differ. Galaxy-wide totals use the first study's 250 pc by 25 pc integral (700,133 clouds) and the class shares of one draw (14.6 / 43.8 / 27.5 / 14.1%).
- Shape fill is the metaball shape's inside test over 6,000 points in the unit sphere for 500 generated clouds; it measures the shape, not whether real clouds are as hollow.
- The real filling of 0.5% to 1% (disc) and 1% to 2% (arms) is from the project's own GEN.47 note and `GMC_ARM_FILLING_FACTOR`, which I did not re-derive; the search did not return a primary number for molecular volume filling in the Milky Way disc.
- The Lynds comparison (about 900 clouds within 1 kpc) uses a layer of 150 pc and a range of 1 kpc that I recalled, not read.
- The O-star count at default scale (270,000) scales the quarter-scale test galaxy's 4,217 by 64 [R]; the mass cut (14 solar masses, GEN.194) lowers it and the hosting ratio makes the nebula count follow.
- Row and byte sizes are estimates from column widths. The 2 KB figure for a stored nebula was not measured on a populated table.
- Screen "ink" counts overlaps and ignores what a frustum culls, so it overstates real coverage; it compares options, not absolute cover.

## Evidence notes

- [S] `galaxy/nebula_field.py`, `tuning.py` (`PHENOMENON_DENSITY_PC3`, `GMC_ARM_FILLING_FACTOR`, `NEBULA_FIELD_*`, `NEBULA_CLASSES`, `NEBULA_HOST_RULES`), `db/models.py` (`nebulae`), `galaxy/nebula_shape.py`, TODO items GEN.152, DB.24, MAP.172 to MAP.182, GEN.202; the cited catalogues and papers in section 2.
- [C] scripts and outputs in `/mnt/project-files/research/nebula-map-visibility/`: `scripts/density_options.py` and `density_options.json`, `shape_fill_fraction.json`, and the first study's `field_counts.json` and `cluster_inner.json`.
- [R] Rice et al. 2016 (1,064 massive GMCs), the secondary Clemens et al. 1991 figure (3.2e5 Bok globules), the 150 pc layer and 1 to 2 kpc range for the Lynds comparison, the 64-times scaling of quarter-scale counts, and the 2 KB per stored nebula.
