# Adaptation traits: how life on a hard world is tougher than life on an easy one

Written the evening of 2026-10-28 EDT by the life chemistry lane (Aux 2, GEN.875). Design only: no code was changed. It feeds GEN.876 (the traits built) and sits inside the life class system ([life-classes.md](life-classes.md), section 7.4): the life of a world is the triple (class, pigment or habitat letter, trait set), and this page is the trait set.

Sources: the life chemistry review (`reports/life-chemistry-review-2026-10-28.md`, section 6 and 6.3, called "the review"), the life class design, and the habitability code on main (`physics/habitability.py`, `habitability_world.py`). Tier labels T1 to T9 and tweak levels L0 to L4 are Boss's ladder (`docs/design/research/stand-in-values.md`). A fact marked (read) was opened in a source for this page or by an earlier lane; (recalled) means from memory or from the review and must be spot-checked before it becomes a hard number in a test. Boss's scope still stands (Boss, 9:12 PM EDT 2026-10-28): open build, if life is possible it develops, most life stops at animals, sentient life is rarest and waits for 9.5.

---

## 1. Short answer

1. A **trait** is one way a body's life has become tougher against one stress that body itself has: radiation, ultraviolet, sulfide, acid or alkali, brine, heat or cold, high carbon dioxide, low oxygen, other toxins. Each trait widens **one named limit** the habitability index already has. Nothing else changes.
2. A trait comes in **three grades**: 1 ordinary (no trait, the limit as it is), 2 tolerant, 3 extreme. Animals and plants are capped below microbes (section 4), because the best known animals are far less tolerant than the best known microbes.
3. A body draws its traits **once, in the life pass**, from its own stress. The more stress, the likelier the trait and the higher the grade; the chance is capped; a body carries at most three. The draws are on the body's own life sub-stream (GEN.872) and are a fixed number of draws, so adding or changing a trait never moves any other draw.
4. **A trait never makes a world habitable to people.** The human scores (PHI-6: Pressure, Temperature, Air, Water, Radiation, Suitability for a captain) never read a trait. A trait widens what **life there** can be, so it moves the life likelihoods (`phi_bio`, `phi_cpx`) and through them the four life scores and the life class scores.
5. A trait is **stored** as a short text (`life_traits`, for example `radiation:2,sulfide:3`), **shown** on the object page ("adapted to: radiation, sulfide") and **read** by the life scores. GEN.876 builds it; it needs the life sub-stream (GEN.872), the biogenic oxygen and timeline work (GEN.855, GEN.750), and the class table (GEN.886), so it goes in the class-switch lock turn as one epoch and one migration with GEN.887, GEN.889 and GEN.890.
6. Three of the old chemistries become traits (life-classes.md 7.4): **Melanin** is the radiation trait (a high dose makes it likelier, as a shield and a repair, not as food), **Zinc-Bacteriochlorophyll** is the acid trait, **Retinal** is a light supplement any class that fixes carbon can carry (section 5.10).

---

## 2. What a trait is not

- Not a species. One body draws one trait set for its whole biosphere (its dominant class); a mixed biosphere is a 9.5 question.
- Not a human tolerance. The human limits in `GAS_LIMITS_KPA`, `PRESSURE_BAND`, `TEMPERATURE_BAND` and the rest stay as they are and keep deciding what a captain can stand.
- Not a ban or a promotion. Where no trait fits the world's stress, the life simply keeps the ordinary limits and the existing scores stand; where the stress is high, a trait only lowers how much the limit costs it.
- Not a new source of energy (except Retinal's small supplement, section 5.10). A tough microbe on a dark irradiated world still needs something to eat; the energy supply term (life-classes.md 7.1 and 8.3) is untouched.

---

## 3. The stress and the limit it moves (the checked table)

One row per stress. "Limit moved" is a name in the code today. "Base" is what the index uses now; the factor is what the grades multiply or add. Every factor and onset is a stand-in with a registry row (section 9).

| # | Stress (the body's own stored value) | Real mechanisms | Limit moved | Grade 2 (tolerant) | Grade 3 (extreme) | Cap for animals and plants |
|---|---|---|---|---|---|---|
| 1 | Radiation: absorbed dose at the surface or habitat, Gy/yr (`dose_gy_yr`) | Fast repair of DNA double-strand breaks and polyploidy (Deinococcus radiodurans); manganese-based antioxidant complexes that protect proteins (Daly and colleagues, recalled); the damage-suppressor protein Dsup of a tardigrade, which halved X-ray DNA damage in human cells that were given it ([R1842], read in a summary); bdelloid rotifers that survive and reproduce after gamma doses about five times those other animals survive (read in a summary; the dose in Gy is recalled as about 1 kGy); melanin shielding in fungi; shielding by water, soil or shade | `DOSE_THRESH_GY` (10 Gy/yr) in `L_rad`, and the separate microbe threshold of section 4.2 | times 10 | times 100 | grade 2 (tolerant) |
| 2 | Ultraviolet: the ground UV index (`uv_index`) | UV-screening pigments (scytonemin, mycosporine-like amino acids), fluorescent down-conversion (the Blue-Optimized idea), living under water or soil (all recalled) | the `UV_BAND` edges used by `uv_score` and the oxygenic penalty on F and A stars | edges times 2 | edges times 4 | grade 2 |
| 3 | Sulfide: the partial pressure of H2S in the air, or the sulfide of the water (`h2s`) | Sulfide-binding haemoglobin and symbiotic sulfide oxidisers in the vent tubeworm; sulfide-quinone oxidoreductase and a sulfide-tolerant cytochrome oxidase in the cave molly Poecilia mexicana, a vertebrate that lives in sulfidic springs (Tobler and colleagues, recalled; the paper was not opened) | the H2S entries of `GAS_LIMITS_KPA` as the life reads them (chronic 0.001 and acute 0.05 kPa for a human) | limits times 10 | limits times 100 (animals) to 1,000 (microbes) | grade 3 at times 100 |
| 4 | Acid or alkali: the ocean pH (`ocean_ph`) | Acidophiles (Zinc-Bacteriochlorophyll's only known host, Acidiphilium rubrum, is an aerobic acidophile, read in the review); proton-pumping membranes; alkaliphiles | the pH window of `L_chem` for the life (the human pH band is untouched) | window wider by 2 pH units each side | wider by 4 | grade 2 |
| 5 | Brine: the water activity of the surface water (`water_activity`) | Halophiles with compatible solutes (retinal pigments are common here); the lowest water activity at which cells divide is 0.605 ([R4222], read) | `WATER_ACTIVITY_CRIT` (0.605) and the chaotropicity limit for the life | animal floor lowered by 0.05 | lowered by 0.10, to the 0.605 record | grade 2; the record is a microbe's |
| 6 | Heat or cold: the surface temperature outside the life's band | Thermophiles and psychrophiles; antifreeze proteins; the warmest growth recorded is 122 C (Methanopyrus kandleri strain 116, read in a summary) | the life's temperature band (the human `TEMPERATURE_BAND` is untouched) | band wider by 15 K each side | wider by 30 K, to the 122 C microbial record | grade 2; animals stay well inside the microbial limits |
| 7 | High carbon dioxide (`pco2`) | Hypercapnia-tolerant physiology of burrowing and deep-sea animals (recalled) | `CO2_TOX_KPA` (2.0, a human number) in `T_metazoa` | times 2 | times 10 | grade 3 at times 5 |
| 8 | Low oxygen (`po2`) | The low-oxygen metabolism of sponges and sessile animals, which live in far less than the 8 kPa band ([R3084], recalled); large surface-to-volume bodies | `METAZOA_O2_LOW_KPA` (8, a human breathing floor, not an animal floor) | floor times 0.25 | floor times 0.05 | grade 3 for sessile grades only (GEN.877 sets the floor by animal grade; this trait is the draw that selects it) |
| 9 | Other toxins: perchlorate, heavy metals, SO2, CO (`so2`, `co`, perchlorate shares) | Metal-resistance genes; perchlorate- and chlorate-reducing microbes | the per-gas limit of that toxin (the human `GAS_LIMITS_KPA` row for it) | limit times 10 | limit times 100 | microbes only (animals take grade 1) |

A trait moves **only** the limit in its row. The table is the list a test can walk: one test per row that a world with that stress and grade 3 scores no lower than the same world with grade 1, and that the human scores of the two are equal to the last digit.

What was checked for this page and what was not: the 122 C record, the halving of X-ray DNA damage by Dsup in human cells, and the relative gamma tolerance of bdelloid rotifers (about five times other animals) were read in summaries; the 0.605 water-activity record and the Zinc-Bacteriochlorophyll host were read earlier by the review. The manganese complexes, the cave molly's enzymes, the grades of the animal caps, the CO2 numbers and the sponge oxygen floor are recalled and are flagged as such in the registry rows. No paper on Boss's cannot-access list is needed.

---

## 4. Grades and caps

### 4.1 Why animals are capped below microbes

The tolerance of the best microbes is orders of magnitude beyond the best animals, and animals have the smallest margins on oxygen, carbon dioxide and heat. So a body whose life reaches animals draws its **animal grade** with the cap of table column "Cap for animals and plants", and a body whose life stops at microbes may draw any grade. The caps are stand-ins (tier T7, tweak L3) set so that a census of generated worlds has some animals on moderately stressed worlds and almost none on extreme ones, which is what the "most life stops at animals" ruling needs.

### 4.2 The microbe radiation threshold is a separate number

Today `phi_cpx = phi_bio * t_metazoa`, so animals and microbes read the same `L_rad` threshold, 10 Gy/yr, which is an ordinary animal's [R0986]. Microbes grow under chronic doses of the order of tens of Gy per hour (recalled; about 5e5 Gy/yr for Deinococcus), so the build adds a **microbe threshold** (`MICROBE_DOSE_THRESH_GY`, proposed 1,000 Gy/yr, T7, L3) that `phi_bio` reads and `phi_cpx` does not, and traits multiply each side's own threshold. That change moves values for every irradiated world, so it ships in the same epoch as the traits.

---

## 5. How a trait is drawn

### 5.1 The stress score

For each of the nine stresses the body has a stress score `s` from 0 to 1, from its stored values:

```
s = clip( log10(stress / onset) / log10(full / onset), 0, 1 )      (the stress rises with its value)
```

For a stress that falls (water activity, temperature on the cold side, oxygen), `stress` is the shortfall against the ordinary limit. The `onset` is the value where the ordinary limit starts to cost the life something; `full` is the value where the extreme grade is needed. Proposed onsets and fulls come from the limit that trait moves (for instance radiation: onset the ordinary threshold, full the extreme grade's) and are registry rows (T7, L3).

### 5.2 Whether the trait is drawn

`p = P_CAP * s`, with `P_CAP` proposed 0.8 (T9, L3): even the harshest world does not always draw the trait, and a world with `s = 0` never does. Two uniform draws per stress, always, whether or not the trait is drawn (a fixed 18 draws per body): `u1` decides presence (`u1 < p`), `u2` decides the grade.

### 5.3 The grade

Every drawn trait is at least tolerant (grade 2). It is extreme (grade 3) if `u2 < GRADE3_SHARE * s^2`, so the share of extreme grades grows with the square of the stress. Proposed `GRADE3_SHARE` 0.5 (T9, L3). The animal caps of section 3 are applied after the draw, so the draw stays the same for any stage.

### 5.4 At most three traits

If more than three are drawn, the three with the highest `s` are kept (ties broken by the table order). Three is a proposed cap (T9, L3); a typical world has none, a stressed one two.

### 5.5 The sub-stream and the order

The draws come from the body's own life sub-stream (GEN.872), after the chemistry or class and pace draws and before the milestone ages, so a world without life draws none and a world with life draws exactly 18 uniforms for traits. A test pins it: changing a trait's table values changes no other draw of any body.

### 5.6 What reads the traits

`habitability.py`'s life functions take a small `Traits` value (a mapping from stress to grade, empty by default). `l_rad`, `uv_score` for the life, the pH window of `l_chem`, the water terms, the life's temperature band, `t_metazoa` (its CO2 and oxygen terms) and the toxin terms multiply or widen the limit they own by the grade's factor. The human `bands()` and all six PHI-6 domains never take a `Traits` argument, which is how rule 4 of the short answer is guaranteed by the signatures, not by care.

### 5.7 Stored, shown, filtered

`life_traits` (TEXT, NULL for none) on planets and moons, in the migration of GEN.876 (with the class columns of GEN.887 and GEN.889), registered in `db/column_units.py` if it is numeric (it is not). The object page lists "adapted to: radiation (tolerant), sulfide (extreme)" under the life section, only when the body has a trait (absent values are hidden). The API carries the text. A map and search filter on a trait is a later, lock-free item.

### 5.8 Order of building

1. The `Traits` value and the limits that read it, with an empty set everywhere: a pure refactor, no value changes (lock-free, can go before the lock turn).
2. The microbe radiation threshold (section 4.2) and the draw (sections 5.1 to 5.5), with the migration: a generation change, one epoch (`--affects fill`), in the class-switch lock turn.
3. The page line, the API field and the filters (lock-free after the migration).

### 5.9 Tests the build needs

- A world with stress and grade 3 scores no lower on every life score than the same world with grade 1, and its human scores equal to the last digit.
- The chance of a trait never falls as the stress rises (monotone) and never exceeds `P_CAP`.
- A world with zero stress draws no trait for any seed.
- The number of draws is fixed: eighteen per body with life, none without.
- No animal draws more than its cap; a microbe-only body may reach grade 3.
- Every trait constant has a stand-in row (the test of GEN.870 extends to the new names).

### 5.10 The three old chemistries

- **Melanin** becomes the radiation trait (the dose is its draw; the byproducts report's "weight times dose factor" survives as the stress score). It is no longer a producer.
- **Zinc-Bacteriochlorophyll** becomes the acid trait. Its only known host is an aerobic acidophile that eats organic carbon, not a photoautotroph (read in the review).
- **Retinal** becomes a **light supplement**: a class that already fixes carbon in a sunlit sea can carry it, and it adds a small share of the light energy (proposed equal to the retinal share of Earth's surface ocean energy budget, a stand-in, T7, L3) without fixing carbon. It is not in the nine-trait table because it moves energy, not a limit (GEN.881).

---

## 6. Tiers and tweak levels

Every number on this page gets a registry row in the same change that builds it, in the style of section 6.18 of `gap-filling-plan.md` (the rows of GEN.870):

| Constant (proposed name) | Value | Tier | Tweak | Why |
|---|---|---|---|---|
| Grade factors per trait (table column 5, 6) | as the table | T7 (shape), T9 (numbers) | L3 | Real mechanisms, planetGen's factors |
| Animal caps | grade 2, grade 3 for CO2, sulfide and low O2 | T7 | L3 | Chosen so the census shows animals on moderately stressed worlds only |
| `MICROBE_DOSE_THRESH_GY` | 1,000 Gy/yr | T7 | L3 | Between ordinary microbes and Deinococcus |
| `P_CAP`, `GRADE3_SHARE` | 0.8, 0.5 | T9 | L3 | Draw shape |
| Max traits | 3 | T9 | L3 | |
| Onset and full per stress | from the limit it moves | T7 | L3 | |
| 122 C record, 0.605 water activity | read values | T2 | L0 | Verified records |

Where a record is a hard physical limit (122 C, 0.605) the extreme grade stops there: the grades widen toward the record but never past it.

---

## 7. What this does to the other items

- **GEN.867 to GEN.881**: unchanged in content. GEN.877 (the animal oxygen floor by grade) and the low-oxygen trait here are the same limit seen from two ends: GEN.877 sets the floor by animal grade, this page's trait decides which grade a body's animals have. Build GEN.877 first.
- **GEN.881** (Retinal and Melanin re-home) is this page's section 5.10.
- **GEN.887 to GEN.890** (the class switch): the traits ride the same epoch and migration as the stored class.
- **TEST.435 / the census** (Bug 2): add the trait counts by stress and grade to the prevalence tables.
- **GEN.855 / GEN.750** (Core 1): nothing here depends on them except that the life pass must exist; the trait draw sits after their draws.

---

## 8. Open questions (not blocking)

1. The animal base for water activity and temperature is stricter than the record the index uses for all life today; the build sets it from the first census (T7, L3) rather than guessing here.
2. Whether a body's trait set should differ between its producers and its animals (a tardigrade next to a lichen) is a 9.5 question; this page gives the biosphere one set.
3. A trait that the body's own stress cannot explain (a shielded cave world on a harsh surface) would need a habitat dose, not a surface dose; the radiation trait reads the dose at the life's habitat (`LETHAL_HABITAT_DOSE_GY_YR` already does this for compact objects), which the build should confirm for under-ice and subsurface life.

---

## 9. Sources and what was not verified

- `reports/life-chemistry-review-2026-10-28.md` sections 6, 6.1 to 6.3 (the table and the proposals), 2.3 and 4.4.
- [life-classes.md](life-classes.md) sections 7.4, 10, 11.2 and 12.
- Read in a summary for this page: [R1842] (Dsup), the 122 C strain 116 record (JAMSTEC release and Wikipedia).
- Read by the review: [R4222] (0.605), [R0986, R1843].
- Recalled, not opened: Daly's manganese complexes, Tobler on the cave molly, [R3084], the CO2 numbers, the sponge oxygen floor. The two pages that failed to load (the sulfide fish paper and the Deinococcus abstract) were not read.
