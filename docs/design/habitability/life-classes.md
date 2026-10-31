# The life class system: what each kind of life eats, drinks, breathes and makes, and how it changes its world

Written the evening of 2026-10-28 EDT, after Boss's 9:10 to 9:13 PM EDT messages (design thread "Design carbon life class system"). Design only: no repo code was changed. It feeds GEN.862 (the design), GEN.863 (the scores) and the sub-items listed in section 12, and it supersedes the scope lines on GEN.874, GEN.875 to GEN.885 where section 11 says so.

Boss's rulings this design follows (all times EDT 2026-10-28):
- 9:10 PM: build in all the life, with what it eats, drinks, breathes and excretes, a class for each kind, "this planet supports Class N life" and "this planet has Class N life at the N stage", so we can plot how each class changes the air. Moved into 8.1 (GEN.862, GEN.863).
- 9:11 PM: "no approvals needed at all, we're building it and a whole new classification system." 9:12 PM: "I don't want to restrict anymore, so open build, the assumption is that if it's possible life will develop." So no approval gates and no chemistry is excluded for being unusual. Only chemistries that cannot exist as a stable liquid on a planet are left out, and section 3.3 says why each one is.
- 9:13 PM: "sentient should be the most rare kind of life with most life stopping at animal life." The stage ladder and the prevalence targets in section 5 make sentient the rarest outcome.
- Earlier, still standing: the share of sentient worlds stays as it is (7:38 PM), GEN.750 means 0.5 / 1.5 / 2.5 / 3 / 60 Gyr, scale 1.0 (9:03 PM), sub-ice mats are not food (GEN.871), Retinal's O2 gate is removed (GEN.861), every constant gets a tier and tweak level (GEN.725, GEN.870), no backward-compat shims, testing mode (new epochs and migrations are fine).

Tier labels T1 to T9 and tweak levels L0 to L4 are Boss's ladder (docs/design/research/stand-in-values.md). Numbers marked (recalled) come from memory or from the two research reports and are not re-opened here; they must be spot-checked before they become a test limit. Source reports used: reports/life-chemistry-review-2026-10-28.md (called "the review"), reports/non-carbon-life-climates-2026-10-28.md ("the climates report"), reports/non-oxygen-life-byproducts-2026-10-28.md ("the byproducts report"), reports/life-scores-2026-10-23.md.

---

## 1. Short answer

1. A **life class** is a way of making a living: what the life drinks (its solvent), what it builds itself from (its base), what it eats (its energy), what it breathes (its oxidant) and what it gives off. Nine classes are defined (section 3); numbers 10 to 19 are kept in reserve, the same way planet class letters are.
2. Every planet and moon gets two statements, kept apart:
   - **"Supports Class N life"**: the classes whose needs the world meets right now, each with a score 0 to 9 (the best one is the headline). It reads only the world's conditions. This is the old habitable index with one more dimension: habitable to what kind of life.
   - **"Has Class N life at the X stage"**: the one class that actually lives there and how far it has got (its dominant life). It reads the star's age and the life timeline. Per Boss's assumption, if a class is supported, that life develops; the stage is how far its clock has run.
3. A **stage ladder** (section 5) is shared by every class, with class-specific words: no life, origin, producers, complex cells, multicellular life, animals in grades (sessile, small mobile, large active), and sentient life last. Sentient is the rarest outcome by design, and most life worlds stop at microbes, plants or animals.
4. Each class has a **gas ledger** (section 4 and 6): what it takes out of the air and puts in, per stage. That is what the stage-aware air (GEN.858) and the biogenic gases (GEN.857) need, and it is what the "Life and air" chart on the planet page plots.
5. **How life is limited is physics, not a ban** (section 7). Each class has a window (temperature, pressure, solvent), an energy supply and a speed. Strange-solvent life (acid, dense CO2, cold methane, ammonia water, silicon in acid) is allowed to develop where its window holds; because its energy and reaction speed are low, most of it ends as microbes or mats and a little reaches simple animals. That is a result of the numbers, which stay tweakable, and not a hard cap.
6. The existing six chemistries are folded in: Chlorophyll a, Blue-Optimized and a new far-red pigment become the three **pigments of Class 1**; Bacteriochlorophylls become Class 2; Retinal becomes a light supplement any class can carry; Melanin becomes a radiation-tolerance trait; Zinc-Bacteriochlorophyll becomes an acid-tolerance trait. The pigment a world gets is chosen by physics (the star's light against the pigment's band), which removes the 15 dead star weights (GEN.867) instead of repairing them.
7. This design sits on Core 1's life timeline (GEN.855 + GEN.750, still waiting for the lock) and Aux 4's life chemistry chain. Section 11 lists exactly what it takes from each and what it changes in their items.
8. Four decisions for Boss are in section 13 (cards in the thread); everything else is a stated design choice he can reverse.

---

## 2. The idea in plain words

A world can be rated for people (PHI-6: how survivable and how pleasant for a human captain). It can also be rated for life of some kind. Today the code does the second one with a single curve that is really "water and carbon life that looks like Earth's". The life class system makes the kind of life a dial.

Each class is defined by six things, in the order a reader wants them:

| Dial | Plain meaning | Stored or derived |
|---|---|---|
| Drinks (solvent) | The liquid its chemistry happens in: water, concentrated sulfuric acid, liquid or dense CO2, liquid methane and ethane, water with ammonia | derived from hydrosphere, temperature, pressure, composition |
| Built from (base) | Carbon for eight classes, silicon-bearing carbon for one | the class |
| Eats (energy) | Starlight in a band, hydrogen from rock, sulfur or iron compounds, radiolytic hydrogen, acid chemistry | derived supply, see 8.3 |
| Breathes (oxidant) | Oxygen, nitrate, sulfate, perchlorate, CO2, acid, or nothing | derived |
| Gives off (waste) | O2, CH4, sulfate, SO2, silica, water and so on | the class ledger |
| Lives where (window) | temperature, pressure, radiation, the air's toxins | stored quantities |

"Supports Class N" is computed from those windows and supplies. "Has Class N life" is the draw plus the timeline.

---

## 3. The classes

### 3.1 The table

Class number is the metabolism. A letter after it is the pigment or habitat variant of that metabolism (the thing the old `LIFE_CHEMICALS` table drew). Tiers are for the claim that the class can exist; tweak levels for the numbers in its window.

| Class | Name | Drinks | Eats | Breathes | Gives off | Needs (window) | Reaches | Tier of the class |
|---|---|---|---|---|---|---|---|---|
| **1** | Oxygenic light-eaters | water | light + CO2 + water, plus N, P, Fe | O2 once it exists (aerobes); before that nothing | **O2** (net, after the sinks fill), buried organic carbon; CH4 only through a small food web | Sunlit water (surface ocean, hycean, or surface water on land); mean habitat 255 to 395 K for microbes, up to about 333 K for animals; CO2 at least about 5 Pa; light in the pigment's band at least 1 percent of Earth's for plants and mats (0.01 percent for deep microbes) (T7, L3) | Animals | T1 |
| 1a | visible pigment (Chlorophyll a) | | 400 to 700 nm | | | star spectrum rich in 400 to 700 nm: F, G, K and bright early M | | T1 |
| 1b | blue and UV tuned (Blue-Optimized, with fluorescence) | | 300 to 500 nm | | | hot bright stars: late B to F; UV shield from fluorescent proteins | | T6 to T7 |
| 1c | far-red (chlorophyll d and f) | | 700 to 760 nm and beyond | | | K and M stars; efficiency penalty for the lower photon energy | | T5 (chemistry), T7 (penalty) |
| **2** | Anoxygenic light-eaters | water | near-infrared light (800 to 1050 nm) + CO2 + a donor: H2, Fe2+ or H2S | nothing (anaerobes); free O2 poisons them | sulfate, sulfur, Fe(III), organic carbon; CH4 comes from their food web | Sunlit water, donor supply, free O2 below about 3 percent (the existing gate); habitat 255 to 370 K (T7, L3) | Microbial mats; animals only if another oxidant source exists | T1 |
| **3** | Hydrogen eaters (chemolithotrophs) | water | hydrogen (geochemical or radiolytic) + CO2 | CO2 (methanogens, acetogens) or sulfate (3b) | **CH4**, acetate, sulfide (3b) | Liquid water in contact with rock that makes H2; no light needed; habitat 255 to 395 K | Microbial | T1 to T2 |
| 3a | geochemical H2 (serpentinisation, vents) | | | | | ultramafic rock plus water plus heat | | T2 |
| 3b | radiolytic H2 and sulfate (deep crust, aquifers) | | | | | radiogenic heat and dose in water-bearing rock | | T2 (Desulforudis, Chivian 2008, read in the review) |
| **4** | Sulfur, iron and nitrogen eaters | water | H2S, Fe2+, NH4+ or CH4 | O2 (4a), nitrate and N-oxides (4b), chlorate and perchlorate (4c) | sulfate, Fe(III), N2, CO2 | Reduced fuel and an oxidant meeting, on vents, seeps, redox layers and under an oxygen sea; oxidant delivery limits the size | Animals (chemosymbiotic tubeworm, clam, mussel and cave-fish style communities) | T2 (4a, 4b), T8 (4c as animals) |
| **5** | Acid-world life | concentrated sulfuric acid (81 to 100 percent) | reduced sulfur, organic carbon, blue and UV light | the acid itself (reduced to SO2 and sulfur); traces of O2 | SO2, sulfur, CO2, water (which dilutes the acid and must be shed) | Very dry world; SO2 and a photochemical route to H2SO4; liquid acid 283 to 373 K (cloud layer) or on warm ground; cloud variant needs a temperate layer under 100 C | Mats and aerial plankton; at most simple colonial or sessile forms | T6 to T7 |
| 5a | cloud-deck plankton | | | | | thick air with a layer at 283 to 373 K and 0.3 to 3 bar | | T6 |
| 5b | pool and crust mats | | | | | warm dry ground, acid pools | | T6 |
| **6** | Dense-CO2 life | liquid or supercritical CO2 | hydrogen (fuel) | CO2 (the solvent is the oxidant) | CH4, water (separates out) | CO2 air or sea: liquid CO2 below 304 K above 5.2 bar, supercritical above that and about 73 bar; very dry | Microbial, at most mats | T7 to T8 |
| **7** | Cold hydrocarbon-sea life | liquid methane and ethane | hydrogen and acetylene (hydrogenation) | nothing | **CH4** | 90 to 190 K; N2 and CH4 air around 1.5 bar; methane at saturation (the existing `solvent()` rule) | Microbial mats | T8 |
| **8** | Ammonia-water life | water with ammonia | hydrogen and chemical gradients in cold alkaline brine | CO2 / carbonate | CH4, N2 | 176 to 273 K water-ammonia oceans (the eutectic is about 176 K, recalled), mostly inside ice moons; ammonia at least 10 percent of the water mass | Microbial | T8 |
| **9** | Silicon-bearing life in acid | concentrated sulfuric acid | as Class 5, with organosilicon chemistry | the acid | SO2, **silica (solid SiO2)** | The Class 5 window plus a silicon-rich crust; silicones resist 98 percent acid at 90 C (read in the climates report) | Microbial | T9 (as a backbone), T7 (as a rare atom) |
| 10 to 19 | reserved | | | | | | | |

Reading the table: "Reaches" is what the physics allows (section 7), not a rule. The classes are ordered by how common the build expects them to be, not by importance.

### 3.2 Which stored values each class reads

All of these exist today except the ones marked NEW, which are the supply quantities of section 8.3.

- Temperature, pressure, hydrosphere (dry, vapour, ice, ice-covered ocean, surface ocean, hycean), `ocean_fraction`, `ocean_ph`, water mass fraction, the ten gas shares (`o2 co2 co n2 ar h2 h2o ch4 h2s so2`), `energy_flux_w_m2`, gravity, radiation dose.
- The star's spectrum in bands (needed for the pigment band times spectrum multiplier, GEN.878).
- NEW: hydrogen supply (serpentinisation and radiolysis), donor supply (H2S, Fe2+), oxidant supply (O2, nitrate, radiolytic oxidants), acid medium (is there liquid sulfuric acid, where), ammonia inventory, liquid CO2 phase. Section 8.3.

### 3.3 Considered and left out, with the reason

Boss's rule is "if it's possible life will develop", so an idea is left out only when the climates report finds that its solvent cannot be a stable liquid on a planet, or the chemistry fails on its own terms. Recorded so no one files them again.

| Idea | Why it is not a class | Source |
|---|---|---|
| Silicon in cryogens (liquid nitrogen, Triton) | The building blocks do not dissolve (silanes at parts per thousand, glycine-sized molecules under 1 ppm); "no cryosolvents are suitable" | Petkowski, Bains and Seager 2020, read in the climates report |
| Sulfur chains, liquid sulfur | Liquid only at 115 to 170 C and the viscosity rises more than 10,000 times above 155 C | climates report 3.1 F |
| H2S, HF, formamide, methanol solvents | Photolysed far faster than they form, or not stable planetary liquids | climates report |
| Pure ammonia seas | Needs ammonia ten times as abundant as water and is photolysed to H2 and N2; the water-ammonia ocean (Class 8) is kept instead | climates report 3.1 D |
| Arsenic in place of phosphorus | The GFAJ-1 claim was refuted | climates report 4 |
| Dusty plasma, nuclear or magnetic matter | No chemistry; fiction | climates report 4 |

Each of these can come back as a new class number if a source changes; the reserve numbers exist for that.

### 10.4 The words (GEN.892)

Built: `src/planetgen/physics/life_class_text.py` with `physics/data/life_class_text.json` gives the name and a sentence for every rung of the ladder (0 to 6, and the animal grades A1 to A3), with words of their own for the rungs each class may reach and generic words for the rest, so a strange-solvent class that does reach a high rung still reads correctly. `lines(class)` and `sentence(class)` give "It drinks ..., eats ..., breathes ... and gives off ..." from the rule file. Nothing shows them yet: the planet page, search and API read them when the class switch (GEN.889) lands.

---

## 4. What each class eats, drinks, breathes and excretes (the ledger)

Per mole of carbon fixed, and the effect on a planet. Gas names use the ten stored gases; "haze" and "colour" are outputs of GEN.857 and the rendering lane. Sizes are the byproducts report's unless marked.

| Class | Takes from the air or water | Puts into it | Net effect on the planet by stage |
|---|---|---|---|
| 1 | CO2 (1 per C), water, light | **O2 (1 per C)**; organic carbon, of which about 0.1 percent is buried and becomes the net O2 source (recalled, T2); CH4 from the food web, falling as O2 rises | Origin: methane and CO2 air. Producers: sinks (reduced volcanic gas, Fe2+, sulfides) soak the first O2; on Earth about 0.3 to 1 Gyr passed before the first rise and about 2 Gyr more before the second (recalled, T2). After the rise: O2 up, CO2 and CH4 down, ozone shield above about 1 percent of today's O2 (recalled), temperature shifts as CH4 and CO2 fall |
| 2 | CO2, a donor (H2, Fe2+, H2S), near-infrared light | organic carbon, sulfate, sulfur, Fe(III); no O2 | Donor drawdown. Food web makes CH4: 100 ppmv to 3.5 percent in anoxic air (9 to 67 Tmol per year), with organic haze once CH4/CO2 reaches about 0.1 to 0.2. Production is donor limited at 1e-4 to 1e-2 of today's biosphere. Never O2. |
| 3 | H2, CO2 (or sulfate for 3b) | CH4, acetate, sulfide | A methane source on open-ocean worlds, small. Sealed sub-ice and aquifer worlds show nothing in the air. Hydrogen supply caps it |
| 4 | H2S, Fe2+, NH4+ or CH4; O2, nitrate or perchlorate | sulfate, Fe(III), N2, CO2 | Draws H2S down, makes sulfate. Size capped by how fast the oxidant arrives |
| 5 | reduced sulfur, organics, the acid | SO2, sulfur, CO2, water | Shifts SO2 against H2S, changes cloud acidity and albedo. On a dry hot world it can drive the SO2 and H2SO4 cycle (the climates report's own bookkeeping, no source gives rates, T6) |
| 6 | H2, CO2 | CH4, water | Methane source; water separates from the fluid (a waste problem) |
| 7 | H2, acetylene (from haze organics) | CH4 | Depletes surface H2 and C2H2, adds CH4. The Titan hydrogen and acetylene anomalies are a lead, not a detection (recalled) |
| 8 | H2, carbonate, chemical gradients | CH4, N2 | Ammonia is destroyed faster than life can use it; the net air effect is small |
| 9 | as Class 5 | SO2, silica | Silica deposits (surface colour, a solid sink); otherwise as Class 5 |

For each class the ledger row has: gas, amount per unit biomass, a stage multiplier (zero before producers, ramping through the producers stage, steady later) and a sink limit. The sinks and the timing for Class 1 are Core 1's GEN.855 numbers (rows LF49 to LF56 and BO1 to BO9 in the registry); this table does not set them. The byproducts report's table is GEN.857's source for Classes 2 and 3.

### Food, water and breath for the captain

Everything above is also what the planet page can say about a world's life, in plain words ("It drinks hydrocarbons and gives off methane"). It sits in one collapsible section that starts collapsed, per the object page rule.

**Is the life edible to people?** Only water-and-carbon classes (1 to 4) count toward the Suitability food rung, as now (Flora is good food, Fauna is food, Boss 10-23). Strange-solvent life (5 to 9) is never food and never counted. Sub-ice and aphotic mats (Classes 3 and 4 microbes) are not Flora and not food (GEN.871, Boss 10-28).

---

## 5. The stage ladder and how rare each rung is

### 5.1 The rungs

The five timeline milestones stay (`tuning.LIFE_MILESTONES`: abiogenesis, photosynthesis, complex cells, multicellularity, technological civilization), because GEN.750 owns them and the database already stores their ages. The ladder gives them class-specific names and adds animal grades from the complexity index (section 7.3), which needs no new timeline step.

| Rung | Milestone behind it | Plain meaning | Class words (examples) |
|---|---|---|---|
| 0 No life | none | conditions are met but life has not started, or conditions are not met | "supports Class N life, none yet" |
| 1 Origin | abiogenesis | chemotroph microbes | all classes: "microbial origin" |
| 2 Producers | the "photosynthesis" milestone, renamed per class | the class's own food maker is running | 1: oxygenating mats; 2: anoxic mats, methane air; 3 and 4: dark chemosynthetic mats; 5: acid cloud plankton; 6, 7, 8: slow microbial communities |
| 3 Complex cells | complex cells | cells with internal structure | protists and plankton |
| 4 Multicellular | multicellularity | macroscopic life | Flora and sponge-grade fauna (Class 1); tube communities (Class 4); colonial mats at best for strange solvents |
| 5 Animals, by grade | multicellularity plus the complexity index | grade A1 sessile or sponge-like; A2 small mobile; A3 large active | grade comes from the energy and oxidant budget (7.3), so a Class 4 vent world can have A1 or A2 and a Class 1 world with high O2 can reach A3 |
| 6 Sentient | technological civilization | species and civilizations | 9.5 (GEN.864): stays as now in 8.1 |

Flora is macroscopic photosynthesis or a surface mat route (Classes 1 and 2); Fauna comes from rung 5 on any class. The rule in `life_scores.scores` that "fauna never exceeds the flora it feeds on" changes to "fauna never exceeds its class's producer base" so chemosymbiotic animals (Class 4) can have Fauna with Flora 0.

### 5.2 Prevalence targets (Boss 9:13 PM EDT)

**A civilization is never forced (Boss, 5:08 PM EDT 2026-10-30).** Adjusting the prevalence, or asking for intelligent life, only tweaks the probabilities: nothing places a civilization, redraws a system until one appears, or squeezes a timeline to make one fit. A request that no draw meets says so plainly.

Sentient is the rarest outcome and most life ends at microbes, plants or animals. The anchor is Core 1's timing review, which Boss chose at 9:03 PM EDT (complex cells 2.5, animals 3, civilization 60 Gyr, scale 1.0): for worlds with oxygenic life, G stars give about 21 in 100 animals and 0.6 civilizations, K about 33 and 1.5, M about 23 and 0.7. Animals to civilizations is therefore about 20 to 40 to 1. These are the bands the census (TEST.435, moved into 8.1) should hold. They are stand-ins, T9 and L3, set from Boss's rule and Core 1's model, not from observation.

| Population | Target band | Why |
|---|---|---|
| Of all worlds that have any life: stop at origin, producers or complex cells | the majority (more than half) | "most life stopping" below sentience; includes every strange-solvent world that cannot get past mats |
| Of worlds with oxygenic life: reach animals (any grade) | about 15 to 35 percent; for Class 1 worlds old enough, about 20 to 35 percent | Core 1's figures, which describe oxygenic life (GEN.907, Boss 6:04 AM EDT 10-29: per oxygenic world, not per life world) |
| Of worlds with oxygenic life: reach sentient | no more than 1 to 2 percent, and always fewer than one in twenty of the worlds that have animals | Boss: rarest; share unchanged from today |

The census (TEST.435) measures the animal and sentient bands on the oxygenic life worlds of systems made around G2V and K5V stars, because the default galaxy holds almost none (4 in 1,246 life worlds on 10-29, chlorophyll chemistries being listed for G and K stars only). The numbers are re-taken and re-baselined after epoch 120 (GEN.867) lands; none are changed before.

Per class, the build thread sets a target band from the class's speed and energy and records it in the registry. Class 1 carries almost all the animals and all the sentient worlds in 8.1. Class 4 animals are chemosymbiotic and small. Classes 5 to 9 are expected to show a low single-digit share of animals at most, and that is an outcome to measure (TEST.435 reports it), not to force.

---

## 6. How a class changes the air, and the chart

### 6.1 What is computed

Core 1's GEN.857 and GEN.858 turn "the kind and stage of life" into a gas change. This design gives them the keys: the ledger rows above are indexed by **class number** (not by the old chemistry text), and the stage by the rungs above. The planet's air at any age is the air before life, plus the sum of ledger rows for its class up to that age, minus the sinks. Class 1 on a standard world follows GEN.855: methane and CO2 at origin, a delay for the sinks, the rise and the ozone shield.

### 6.2 The "Life and air" chart on the planet page

A small line chart under the Life section: x is age (Gyr) from 0 to the system's age; one line per major gas (O2, CO2, CH4, plus H2S or SO2 where the class moves them); stage boundaries are marked as ticks (origin, producers, complex cells, multicellular, animals). It is computed on demand from the stored milestone ages and the class ledger, so it needs no stored column. It answers Boss's "plot how each would change their environment". A second view, for the census and the galaxy map, compares classes: median O2, CH4 and haze of Class 1, 2 and 3 worlds at a stage (a table in the census report, not a page, in 8.1).

### 6.3 What the plot does not claim

Rates for Classes 5 to 9 are the climates report's redox bookkeeping, tier T6 to T8, L3 to L4. The chart marks them "speculative chemistry" with the tier, as the climates report recommended.

### 6.4 As built (GEN.891, GEN.857, GEN.858, GEN.859)

The gas ledger of section 4 is `physics/data/life_ledger.json`, read by `physics/life_air.py` and applied in the life pass after the oxygen (`generation/life.py`, `_build_air`):

- **Keyed by class number.** One row per class (1 to 9) with a status. All nine are `built` (GEN.894 built 4, 5, 7 and 9). Every number and row carries its tier and tweak level, with summary rows LF60 and, for the sulfur and acid rows, the same file.
- **Methane.** The flux is the class's fuel supply (Class 1: the organic carbon its photosynthesis buries, half a methane each; Class 2: `donor_supply` times Earth's 10 Tmol a year; Classes 3, 6, 8: `h2_supply` times Earth's 2.3 Tmol a year of hydrogen over 4 H2 per CH4) over the world's area. The air holds the flux times the 30,000 year life of methane in an oxygen-free air, at most 3.5 percent of its molecules, and the methane's carbon comes out of the air's own carbon dioxide (10 percent of the carbon stays carbon dioxide), so no world makes methane from carbon it does not have. A gain under 1 Pa (`biogenic_oxygen.MINIMUM_PA`) leaves the air as it was.
- **Hydrogen.** The classes that eat hydrogen (2, 3, 6, 7, 8) draw it down to 1e-4 of the molecules; Class 7 turns each mole it takes into two thirds of a mole of methane (C2H2 + 3 H2 -> 2 CH4; the acetylene itself is not an air gas).
- **Hydrogen sulfide (GEN.894).** The sulfur eaters (4) and the acid eaters (5, 9) draw hydrogen sulfide down to 1e-8 of the molecules (a design floor, T9 L4); Class 4 makes sulfate in the water, not a gas, while Classes 5 and 9 give back the same moles as sulfur dioxide, so the air shifts from H2S to SO2. The caps are taken before the products, so a gas never ends above its cap.
- **Stage.** Nothing before the second step of the schedule (the producers) and, for Class 1, nothing once the oxygen has risen or while the air holds free oxygen (more than 0.001 of its molecules): the oxygen has burned the methane (`biogenic_oxygen.burn`). A world with no life keeps its abiotic air.
- **Redone from the new air.** The haze (`aerosols`), the ultraviolet shield and dose (`radiation`), the colour and the habitability score are generated again from the air the life left. The stored surface temperature follows the new greenhouse (`biogenic_oxygen.store_air`).
- **Not modelled.** The cooling of the ground by the haze, the hydrosphere of a world the new air moves across the freezing point (its stored state of water stays, as after the oxygen), the sulfate, sulfur and silica the sulfur and acid classes leave in the water and rock, and the acetylene of Class 7.

---

## 7. How life is limited without a ban

### 7.1 The three limits every class has

1. **Window**: temperature, pressure, solvent stability, radiation and toxins (section 3 and `L_*` terms in the habitability index).
2. **Energy supply**: how much usable energy reaches the biosphere. Light classes use the star's flux in the pigment's band (GEN.878). Chemical classes use the supply of fuel and oxidant (section 8.3). A class's energy supply sets its biomass, and its biomass sets how fast it fills a surface and how loud its air signature is.
3. **Speed**: reaction rate against temperature. The climates report's plain rule is a Q10 of 2 to 3 per 10 degrees (Petkowski 2020, read there), with diffusion-limited steps in cryogenic liquids staying fast (relative rate 0.35 to 0.91, Bains 2004, read there). Core 1's GEN.750 clock does not use that rule: it is a Boltzmann factor (E = 0.30 eV, reference 288 K) clamped to 0.25 to 2.5 (read in Core 1's design doc), so a cold world is at most four times slower, not 1e8 times. This design keeps Core 1's clock for the first rungs and lets the energy supply and the oxidant (7.3) decide the upper rungs.

### 7.2 Why strange-solvent life ends up small

The aerobic yield is about nine times the anaerobic yield (Gale and Wandel, read in the review). Acid, CO2 and cold-methane biospheres have far lower energy density than oxygen respiration, and Titan's chemistry would run about 1e8 times slower under a plain Q10 rule (the climates report's arithmetic). Under Core 1's clamped clock the slowdown is at most four times, so the speed term alone does not stop a cold class; its small energy supply and lack of an oxidant do. Their upper rungs therefore mostly never arrive. That is the physics doing the ceiling; the tuning knobs stay open (section 13, card 1 and 2).

### 7.3 The complexity index and animal grades

A class reaches a rung when its **complexity index C** passes a threshold. C is the class's biomass energy (relative to Earth's biosphere = 1) times the oxidant yield factor (1 for aerobic O2 respiration, about 0.11 for anaerobic) times the temperature-speed factor. Starting thresholds (T9, L3; fitted on a first probe of the chemical classes and to be reset from the census): complex cells 1e-5, multicellular bodies 3e-5, sessile or sponge-grade animals 1e-4, small mobile animals 1e-3, large active animals 1e-1. This is the generalisation of GEN.877 (the animal oxygen floor by grade): for Class 1 the oxidant factor reads the world's own O2, with the sponge floor of about 0.1 to 0.8 kPa from Mills et al. 2014 (recalled), rising with body size by the diffusion argument of Catling, Glein, Zahnle and McKay 2005 (recalled). It also gives routes R2 to R5 a home: R2 is Class 4a, R3 is a low-O2 Class 1 world at grade A1, R4 is a rare anoxic animal niche on Class 2 and 3 sulfidic worlds, R5 is Classes 4b and 4c. They are no longer separate items with an approval gate.

### 7.4 Traits stay, and move thresholds

GEN.875 and GEN.876 (adaptation traits: radiation, UV, sulfide, acid, brine, heat or cold, CO2, low O2, toxins) are unchanged in kind. A trait widens a class's window, never makes a world habitable to people. Retinal (a light supplement for any class that has a carbon fixer) and Melanin (radiation tolerance, dose in its draw, GEN.860 and GEN.881) join the trait list; Zinc-Bacteriochlorophyll becomes the acid tolerance trait. So the **life class is the triple (class, pigment or habitat letter, trait set)**, exactly as the review's section 6.2 predicted.

---

## 8. Scoring, supply and dominance

### 8.1 The score of a class

For each class c, the score S_c from 0 to 9 is `ceil(9 * L_c)` with the same shape as today's PHI_bio: a Liebig blend (minimum weight 0.65, `LIEBIG_MIN_WEIGHT`) of

- `L_solv`: the class's solvent exists and is stable (liquid water, acid, CO2, methane...);
- `L_temp`: the habitat inside the class's temperature window (smooth falloff outside, as now);
- `L_press`: pressure where the class needs it (a floor for liquid water, a range for dense CO2, a cloud-layer rule for 5a);
- `L_ener`: the class's own energy supply on the log curve `ENERGY_FLUX_MIN` to `ENERGY_FLUX_REF`, with light classes reading the pigment's band times the star's spectrum (GEN.878) and chemical classes reading their fuel supply;
- `L_redox`: oxidant and donor present for the classes that need them (4, 2, 3b, 5), and free O2 below the class's toxic limit (2, 3);
- `L_rad` and `L_chem`: dose and toxins against the class's own limits (microbe limits for microbial rungs, animal limits for animal rungs), with traits widening them.

**Class 1a on an Earth-like world is exactly today's PHI_bio**, so the existing habitable index becomes the near-Earth carbon and water case with no change to its values (Boss, 10-28: "the existing single habitable index stays as the near-Earth carbon and water class").

Complex-life score per class: `S_cpx_c = S_c * T_c(index C)`, the generalised `t_metazoa` of 7.3.

### 8.2 Dominant class (the "has" side)

Every class that scores at least `HABITABILITY_LIFE_SCORE_MIN` (0.005) develops (Boss: if it is possible, it will develop). The one reported as "has Class N life" is the class with the largest **net primary production** among the supported classes. Net primary production is the energy supply times the class efficiency. Where oxygenic and anoxygenic life both fit a lit ocean, Class 1 wins by orders of magnitude, as on Earth; Class 3 wins under ice because nothing else fits; Class 5 wins on a dry acid world. The other supported classes are listed as "also supports" and keep their scores. The class that lives on a world is always among the classes stored as supported: when its own life has used up what it lived on (the donor or hydrogen it drew down, the air it changed) and the conditions after the air build score it 0, it is stored at the score it was chosen at (TEST.486). A tie within a factor of 2 is broken by a weighted draw (one of the few draws in this design; the pigment letter within Class 1 uses the same rule).

This replaces the star-letter weights in `LIFE_CHEMICALS` and `STAR_EVOLUTION` (the 15 dead weights of GEN.867). The old weights were T4, L3 guesses; the physics rule is stronger, and it makes the M-dwarf far-red route (1c) automatic.

### 8.3 New supply quantities

New derived and stored numbers are the main build cost of the chemical classes. Each needs a units registry entry and a stand-in row.

| Quantity | Meaning | Feeds | Note |
|---|---|---|---|
| `h2_supply` | hydrogen from serpentinisation and from radiolysis | Classes 3, 6, 7, 8 | from rock type (GEN.441 crust composition), water, interior heat and radiogenic heat; tier T6 |
| `donor_supply` | H2S and Fe2+ delivered | Classes 2, 4 | from volcanism, ocean chemistry and the existing H2S and iron stand-ins |
| `oxidant_supply` | O2, nitrate, radiolytic oxidants, perchlorate | Class 4 and animals on aphotic worlds | the cap of GEN.880; no Europa O2 flux number from second-hand sources |
| `acid_medium` | liquid sulfuric acid present: 0 none, 1 a cloud layer, 2 a pool | Classes 5, 9 | built (GEN.888): `physics/acid_medium.py`, the supply of sulfur dioxide, an oxidiser and water, the water activity against the acid strength (81 percent and up) and the liquid window of 283 to 373 K; Class 5a reads 1 or more, 5b and 9 a pool (2). Registry row LF61 |
| `sulfur_supply` | reduced sulfur for the acid eaters | Classes 5, 9 | the air's sulfur dioxide over Venus's 150 ppm, at most 1 (GEN.894, `physics/life_supply.py`); needs no water in contact with rock |
| `silicon_rich_crust` | the rock's molar Mg/Si at or under 1 | Class 9 | `physics/strange_solvents.py` (GEN.894), registry row LF62 |
| `co2_phase` | liquid or supercritical CO2 at the surface | Class 6 | from temperature and pressure against the CO2 phase diagram |
| `ammonia_water_mass_fraction` | ammonia as a share of the water's mass | Class 8 | built (GEN.894): a tenth of the non-water ice over the water (`physics/strange_solvents.py`, registry row LF62); the weakest number of the class |

---

## 9. How it ties to the habitability index and the four life scores

- **PHI-6 stays**, the human rating. Its Habitability dimension (`phi_bio`, `phi_cpx` inside it) becomes "the best score over the water-carbon classes 1 to 4 for the life that is on the world", so nothing about PHI-6 moves for worlds with Class 1 life.
- **The four life scores (Microbial, Flora, Fauna, Sentient)** keep their meaning and their 0 to 9 scale (GEN.649). They are read from the world's **actual** class and stage. Microbial from rung 1, Flora from a photic or surface-mat producer (Classes 1 and 2), Fauna from the animal grades of any class, Sentient unchanged in 8.1. Flora under ice and in aquifers is 0 (GEN.871).
- **The Suitability food rung** (`life_scores.food_stage`) counts only Classes 1 to 4 life, and not Class 3 or 4 microbes under ice.
- **A planet has up to nine class scores** (one per class), plus the headline "supports Class N" (the best), "also supports" (the others at or above the line) and "has Class N at the X stage".

---

## 10. What gets stored, shown and filtered (GEN.863)

### 10.1 Storage

Planets and moons each gain:
- `life_class` (the actual class, small integer, NULL for none) and `life_class_variant` (a short text like `1c`), replacing `life_chemical` (no compatibility shim: the old column is dropped in the same migration, as `phi_bio` and `phi_cpx` were in migration 0136; Boss's rule). The population, search and map readers move to the new column in the same change;
- `life_class_best` (the headline supported class) and `life_class_best_score`;
- one table `life_class_scores(object_id, life_class, score)` for the per-class scores (nine rows at most, only those at or above the line), indexed on class and score so the map filter and the search facet are cheap. The migration also stores units in `db/column_units.py`.

`life_stage` stays as the milestone key; the rung and animal grade are derived from it and from the world's complexity index at read time (the grade is also stored as `life_animal_grade` TINYINT if reading it proves slow; the lane decides).

### 10.2 Display

- Planet and moon page: a Life section under the at-a-glance block. Lines, in this order: one paragraph (what lives there with its stage and animal grade, each kind of life with a prevalence word, the classes the conditions suit with their suitability), then a collapsible block "What it eats, drinks, breathes and gives off" with the class's lines from the ledger, then the Life and air chart, then the traits ("adapted to: radiation, sulfide"). Absent or zero values are hidden unless absence matters (the page rule).
- API: the same fields as numbers and codes; a small API version bump (testing mode allows it); the API compatibility page (API.4) lists it.
- Map and search: a world filter "supports Class N" and "has Class N" with a minimum score, using the table index; the sector colouring by habitable locations (MAP.143) is unchanged.
- System page and lists: the class number beside the four life icons.

### 10.3 Rule file and code

Mirroring the planet class classifier: a rule file `src/planetgen/physics/data/life_classes.json` and a pure module `physics/life_classes.py` with `classify_life(values)` returning the per-class scores, built on the same `holds` and `quantities` machinery as `physics/planet_classes.py`. A shadow sweep (like `class_shadow_cli.py`) can run the classifier on generated bodies with no stored change before the draw switches over. Every window number lives in the file with its tier and tweak level; the docs page is generated from the file, as the planet class pages are.

---

### 10.5 As built (generator epoch 120, schema 175; GEN.867 to GEN.872, GEN.887, GEN.889, GEN.890)

The draw switched over in one lock turn. What differs from the plan above:

- **Storage.** `life_class` is the class number (TINYINT) with `life_class_variant` (CHAR(1)), `life_class_best`, `life_class_best_score`, `life_animal_grade` and `life_complexity` are stored on the body, and the scores are two tables, `planet_life_class_scores` and `moon_life_class_scores` (a child table for each body kind, like the gas tables). `life_chemical` is gone and `species.life_chemical` became `species.life_class` (the id). Migration 0175 maps the six old chemistries to a class and works nothing else out; the scores, grade and complexity come when a galaxy is brought to epoch 120.
- **The supplies** (section 8.3) are not stored: `physics/life_supply.py` works the six shares of Earth's (`h2_supply`, `radiolytic_h2_supply`, `donor_supply`, `nitrate_supply`, `perchlorate_supply`, `hydrogen_acetylene_supply`) from the columns the body stores, with every number in `physics/data/life_supply.json` carrying its tier and tweak. `nh3` and the acid medium are still GEN.888 and GEN.894.
- **Dominance** (8.2) is the supported class with the largest net primary production (`tuning.LIFE_NPP`); classes within a factor of 2 are tied and an NPP-weighted draw picks among them. The draw is made on a stream of the body's own, seeded by one draw of the run's (GEN.872), so a change to this model no longer moves the names and flavour text of later bodies.
- **Complexity** (7.3) is `physics/life_complexity.py`: C is the class's net primary production times the yield of its oxidant times the temperature rate factor, the thresholds, the oxygen floors for an oxygen breather and the anoxic route R4 are `tuning.LIFE_COMPLEXITY`, and the highest milestone a body can reach is capped by C. Flora needs a class that eats light (GEN.871); Fauna needs an animal grade and never exceeds its producer base.
- **Text** (GEN.869, GEN.892): the life and timeline sentences name the class ("Class 1a life, the oxygenic light-eaters ...") and use the class's own words for each rung and for the animal grade.

### 10.6 As built: the planet page's Life section (UX.237, first part)

The Life section of a planet's or moon's page (`web/item_pages.py`, `web/lib/planet_layout.life_class_view`, `item.html`) reads the system-detail row (API.37: `life_class_id`, `life_class_scores`, `life_class_best`, `life_class_best_score`, `life_stage`, `life_animal_grade`) and, for a world with a living class, the life readings above (`GET /api/planets/<id>/life-readings`). At a glance, under the table, is one Life paragraph (`life_class_view(...)["glance"]`): the sentence of each kind of life that is there ("Microbes: abundant, stable, secure. Flora: occasional, declining, precarious."), or "Nothing lives here yet." when the conditions suit a class and none has appeared, then "Conditions suit Class 1a (oxygenic light-eaters) best, then Class 2 (...) and Class 3a (...)." (the best-suited class and at most two others, then "others"; no scores). A world that could not have life shows no paragraph. The summary no longer lists life scores. Under the "Life" heading the section opens with the lead paragraph (`["paragraph"]`): "Life here is Class 1a, at the stage of small mobile animals (animal grade A2)." (the class that lives there, with the stage names of `life_class_text` and the animal grade only at the multicellular stage; "Nothing lives here yet." when the conditions suit a class and none has appeared) and "Conditions here suit Class 1a life (oxygenic light-eaters) best, suitability 7 of 9; they also suit Class 2 (...), Class 3a (...)." (the best-suited class, then the others, best first, with the stored scores). Then the table of the kinds (`["readings"]`: kind, how much with the score out of 9, trend, security and, when known, how long the conditions last), then the stored facts (`["facts"]`: stage, animal grade, pace of evolution, age of the system and the age the stage was reached), then what the life looks like (the description's `life` prose: its colors in and outside visible light), then the collapsed block "What it eats, drinks, breathes and gives off" (the class's `drinks`, `eats`, `breathes` and `gives_off` from the rule file; for a world that supports a class but has none, the best-supported one). The description keeps the sentences "Class 1a life, the oxygenic light-eaters, lives on this planet, evolving at a standard pace." (`life-class`), "The system is about 4.6 Gyr old." (`timeline`) and "Its latest milestone is ..." (`milestone`) as prose blocks for the text and Markdown forms; the page leaves them out because the paragraph, the table and the facts say them. The stage sentences name the whole class ("Class 3a life has grown ...", `life_class_text.stage`). A line with nothing to say is left out, and a body that supports no class and has no life shows none of them. The class number (`1a`) follows the four life icons in every row that has the icons and a `life_class_id` (`habitability_symbols.life_row`: the object page and the system page). A class stored on a body with no life stage (its star is younger than the first step of its schedule, so the world suits the class but no life has started) is not the class that lives there: `life_class_id`, the class markers on the Galaxy Map, the `life` filter and its column, the dominant species line and the system map side panel all leave it out (GEN.1001); the description says the conditions suit the class and none has appeared yet. The Life and air chart (UX.237, second part) follows the words block, collapsed: `GET /api/planets/<id>/life-air` (and `moons`) gives `physics.life_air.air_history`, which works the air backwards from the stored one and stores nothing. Before the producers (the second step of the schedule) it is the air as the life found it (the methane's carbon back in the carbon dioxide, a class that splits water without its oxygen); from the producers on it is the air with the class's methane and no oxygen (the ledger's own `methane_flux_mol_yr` and `methane_target_mol`); from the oxygen rise on (`biogenic_oxygen.Profile.delay_gyr` after the producers) it is the air now. A class whose methane does not end at the rise keeps one air from the producers on. The page draws oxygen, carbon dioxide and methane as partial pressures on a log axis (`web/lib/life_air_chart.py`, an SVG with the stages numbered, the same figures in a table under it); a class whose ledger tweak is L3 or L4 is marked "Speculative chemistry" with its tier and tweak. Not drawn: the sulfur gases, the hydrogen the hydrogen-eaters draw down, the reduced gases the rise burned (an abiotic air's reductants are not stored), and a world the life did not change. Not yet: the adaptation traits (GEN.875 and GEN.876 store no traits yet).

## 11. Overlap with Core 1 and Aux 4, and what changes in their items

Core 1's life timeline design (branch version, unmerged) was copied to /mnt/project-files/design-drafts/life-timeline-and-oxygen-2026-10-28.md at 9:19 PM EDT and I read it, with Core 1's answers to my three questions. GEN.855 + GEN.750 wait for the lock. What I took from it: five hard steps for every chemistry (means 0.5, 1.5, 2.5, 3, 60 Gyr, scale 1.0); a temperature clock clamped 0.25 to 2.5; a first-light scale by metabolism (`LIFE_FIRST_LIGHT_SCALE`: retinal 0.4, anoxygenic 0.6, radiotrophic 1.0, oxygenic 1.0); oxygen made only by oxygenic life in a sunlit sea (`SUNLIT_OCEAN`: surface or hycean ocean), none under ice; the only gas ledger today is oxygen, keyed by `evolution.releases_oxygen`.

**The one real conflict.** In Core 1's timeline, complex cells need free oxygen, so a world whose life is not oxygenic and whose air has no oxygen stops at the producers rung. Under this design that gate must read the **class's oxidant** instead of free O2: oxygen for Class 1 and 4a, nitrate or perchlorate for 4b and 4c, sulfate or CO2 for the anaerobic classes at the low anaerobic yield (about 0.11 of aerobic), so Classes 2 to 9 can pass the producers rung slowly, as Boss's open build and the review's routes R2 to R4 need. This is a change to GEN.750's third step, small in code (swap the oxygen test for a class lookup of the oxidant), and it belongs in L5 or L6 in section 12, after GEN.855/750 merge, not in their change.

### 11.1 With Core 1 (owns GEN.855/750 and GEN.857 to GEN.861, 872, 873)

- GEN.855 + GEN.750 merge first, as planned. Nothing in this design changes them. Core 1 confirmed the first-light factor is keyed by the `metabolism` field and can be re-keyed by class by swapping the lookup.
- GEN.857 and GEN.858 (ledger and stage-aware air): ask Core 1 to key the ledger rows by **life class number** and the stage by the rung, as in section 6.1, instead of by chemistry text. This avoids writing the ledger twice. Core 1 confirmed the only ledger today is oxygen, keyed by `releases_oxygen`, and that re-keying GEN.857/858 rows by class is a swap of the lookup in two places.
- GEN.861 (chemoautotroph carbon base for Retinal and Melanin worlds) is Class 3 and 4 in this design: a chemotroph base with Retinal as a supplement trait. It is consistent with Boss's Retinal ruling; the content is the same, the home changes.
- GEN.872 (the per-body life sub-stream) is required before the draws of Section 8.2 change. Core 1's call whether it rides the GEN.855 epoch.
- Epochs: the draw switch (8.2) and the new columns (10.1) are one epoch and one migration. The classifier alone (10.3 shadow) is lock-free.

### 11.2 With Aux 4 (owns the life chemistry chain GEN.867 to GEN.884, TEST.436)

| Item | What this design does to it |
|---|---|
| GEN.867 (star lists agree with weights; no chemistry on a non-water solvent) | Replaced: the weights go away (8.2). The "never a water chemistry on a non-water solvent" guard becomes the class's solvent term. Re-scope to "remove the star weights and star lists; add the solvent term" |
| GEN.868 (chemistry only after abiogenesis) | Kept as is; it becomes "class only after abiogenesis" |
| GEN.869 (text and labels) | Kept; the text moves from chemistry text to class text (stage names, eats, drinks, breathes, gives off) |
| GEN.870 (registry rows for every life constant) | Extended to every class window and the complexity thresholds; the lock-free part (rows for today's constants) goes ahead now |
| GEN.871 (sub-ice mats scored as Flora) | Kept; Class 3 and 4 microbes under ice are not Flora (5.1) |
| GEN.873 (standing tests) | Extended: classes in place of chemistries |
| GEN.874 (non-carbon proposal) | **Retire**: superseded by Classes 5 to 9 and the "left out" list (Boss 9:11 PM EDT) |
| GEN.875, GEN.876 (adaptation traits) | Kept; the trait list gains Retinal, Melanin and the acid trait (7.4) |
| GEN.877 (animal O2 floor by grade) | Generalised to the complexity index (7.3) |
| GEN.878 (chemistry-specific energy) | Kept; it is the pigment band times the spectrum, now the pigment selector of Class 1 and 2 |
| GEN.879 (far-red oxygenic) | Kept; it is Class 1c |
| GEN.880 (chemolithoautotroph producers) | Kept; Classes 3 and 4, with the supply quantities of 8.3 |
| GEN.881 (Retinal and Melanin re-home) | Kept; section 7.4 |
| GEN.882, 883, 884 (routes R2 to R4) | Kept as the animal routes of 7.3 (R2 is 4a, R3 and R4 are Class 1 and 2 low-oxidant cases) |
| GEN.885 (R5 proposal) | **Retire**: Classes 4b and 4c, no approval needed |
| TEST.436 and TEST.430 | The class census (TEST.435, moved to 8.1) lives in Bug 2's suite |

Aux 4's chain order (GEN.870, 869, 875, then the lock bundle 877, 878, 879, 880, 876, 881, 882 to 884) already puts the pieces this design needs in a workable order. The one dependency to add is that Aux 4's Class 1 variants (879) and chemolithotrophs (880) should be written as class-table rows from the start (section 10.3), not as new entries in the old `LIFE_CHEMICALS` dictionary; otherwise they will be rewritten.

### 11.3 Sequencing and lock turns

1. Now, lock-free: this design merged as a docs page (DOC), GEN.870 rows, the shadow classifier and rule file (A in section 12), the text work.
2. After GEN.855/750: Core 1's GEN.857 to GEN.861 on the class-keyed ledger.
3. One lock turn for the class switch (stored classes, scores, draw replacement, migration), carrying GEN.872 and Aux 4's bundle where they fit.
4. A second lock turn for the supply-heavy classes (3 to 9) when their supply quantities exist.
5. Display, API, map filter and the Life and air chart are lock-free after the migration.

Suggested owners (TODO decides): the classifier, rule file, text and display to Aux 4's lane because it already holds the life chemistry chain and does not retire; the ledger and air to Core 1; the census to Bug 2. Note Boss's standing rule that no new Aux or Math work starts, so if Aux 4 retires before this is built, the class work goes to Core 1 after GEN.861.

---

## 12. Sub-items for TODO to file

Numbers are TODO's. "Needs" are prerequisites. All are 8.1.

| # | Item | Needs | Lock? |
|---|---|---|---|
| L1 | GEN.862 re-scoped: "Design: the life class system, merged as docs/design/habitability/life-classes.md from the design draft" (this note), and change its text to say Boss ruled 9:10 to 9:13 PM EDT that all researched chemistries are in and sentient is the rarest rung | none | no |
| L2 | The life class rule file and classifier: `physics/data/life_classes.json`, `physics/life_classes.py`, the shadow sweep, the generated docs page | L1 | no |
| L3 | Supply quantities: `h2_supply`, `donor_supply`, `oxidant_supply`, `co2_phase`, `nh3_inventory` and their stand-in rows, with units registry entries (builds on GEN.880 and GEN.441) | L1 | yes (new stored values) |
| L4 | The acid medium model: a sulfuric cloud layer and acid pools from the temperature profile, SO2 and water shares (nothing models it today); fills the dead `solvent()` branch | L1, L3 | yes |
| L5 | Class switch: `life_class`, `life_class_variant`, `life_class_best`, the `life_class_scores` table and migration, drop `life_chemical` and the star weights and lists, dominant class by net primary production (GEN.863 build part 1) | L2, GEN.872, GEN.855, GEN.750 | yes |
| L6 | Complexity index and animal grades per class; the complex-cells gate reads the class's oxidant instead of free O2 (a change to GEN.750's third step); the Fauna rule for chemosymbiotic life (generalises GEN.877) | L5, GEN.877 | yes |
| L7 | The class ledger and class-keyed stage-aware air (tell Core 1; shares GEN.857 and GEN.858), including Classes 2 to 9 | L5, GEN.857, GEN.858 | yes (in Core 1's second change) |
| L8 | Strange-solvent classes 5 to 9: windows, supplies, texts, stand-ins and tiers | L3, L4, L5 | yes |
| L9 | Planet page: Life section, eats-drinks-breathes-gives-off block, Life and air chart, class number on lists | L5, L7 | no |
| L10 | API fields and version bump; map and search filters for "supports Class N" and "has Class N" (GEN.863 build part 2) | L5 | no |
| L11 | Class text: stage names and descriptions per class replacing the chemistry text (GEN.869, GEN.881) | L1 | no |
| L12 | Census: prevalence per class and stage inside researched ranges, sentient the rarest rung, animals to civilizations at least 20 to 1 (TEST.435 moved from 9.5 to 8.1; Bug 2 owns, extends TEST.430) | L5 | no |
| L13 | Registry rows for every class window, supply and complexity threshold with tier and tweak level (extends GEN.870) | L2 | no |
| L14 | Retire GEN.874 and GEN.885; re-scope GEN.867, GEN.873, GEN.877, GEN.879, GEN.880, GEN.881 as in 11.2 | L1 | no |

GEN.863 itself ("Life class scores... shown on the object page, the API and the map filters") is L5 + L9 + L10. GEN.864 and GEN.865 stay in 9.5.

---

## 13. Decisions for Boss (cards)

Answered (EDT 2026-10-28): card 1, option 1, physics decides with no hard cap (9:19 PM); card 2, option 1, origin uses the faster rate (9:19 PM; note Core 1's clamped clock, 7.1, already makes this the default); card 3, option 1, most food production wins (9:20 PM). card 4, option 1, Classes 1 to 4 first and 5 to 9 second, both in 8.1 (9:20 PM). All four cards are answered and the design stands as written; sub-items L1 to L14 (section 12) are final, with L8 (Classes 5 to 9) as the second build.

1. **How far can strange-solvent life go (acid, dense CO2, cold methane, ammonia water, silicon in acid)?** Option 1 (recommended): physics decides. Its energy supply and reaction speed set the ceiling, with no hard cap, so most of it ends as microbes or mats and a few percent reach simple animals. Option 2: hard cap at colonial or sessile animals. Option 3: same ladder as water life.
2. **Cold methane worlds are about 100 million times slower by the plain temperature rule.** Should life start there at all? Option 1 (recommended): yes, the origin step uses the faster diffusion-limited rate (Bains 2004), later steps stay slow, so they end as microbes. Option 2: plain temperature rule, so they "support" life but stay empty.
3. **When several kinds of life would fit a world, which is the dominant one?** Option 1 (recommended): the one with the most food production wins and the others are listed as "also supports". Option 2: a weighted random pick.
4. **How much of the system goes into 8.1?** Option 1 (recommended): the store, the display and Classes 1 to 4 first, then the strange-solvent classes 5 to 9 in a second build, both inside 8.1. Option 2: all nine classes in one go. Option 3: only Classes 1 to 4 in 8.1, the rest in 9.5.

One more choice I made without a card, listed so it can be reversed: the pigment a world gets (Class 1a, 1b, 1c) is picked by physics (star light against the pigment's band) instead of the old star-letter weights.

---

## 14. Sources and what was not verified

Read for this note: the three reports named at the top; docs/design/habitability/habitability-index.md section 5 and 12, condition-based-planet-classes.md; `tuning.py` lines 2174 to 2420 (LIFE_CHEMICALS, STAR_EVOLUTION, EVOLUTIONARY_TIMELINES, LIFE_METABOLISMS); `physics/life_scores.py`; `physics/data/planet_classes.json` and the quantity list of `physics/planet_classes.py`; `db/models.py` life columns; docs/TODO.md items GEN.855, 750, 857 to 885; docs/plan/releases/version-9-5-life-classes.md; Core 1's and Aux 4's recent thread status; main at b851509d8.

Core 1's design doc was read (the branch copy of 9:19 PM EDT); its code was not. No probe was run for this design (the review's measurements are quoted, not repeated).

Recalled, not re-opened: the Class 1 timing for the first and second O2 rises, the 0.1 percent organic burial, the ozone shield threshold, the water-ammonia eutectic near 176 K, CO2 critical constants (304 K, 73 bar), the sponge O2 floor (Mills 2014), the Catling et al. 2005 size argument, the CO2 compensation of about 5 Pa and the light compensation thresholds. The proposed complexity thresholds in 7.3 are starting values, not sourced numbers. Spot-check every number before it becomes a test limit, as the review asked.

## Life readings: abundance, trend and security (GEN.1210)

`physics/life_readings.py` reads three words for each kind of life (microbes, flora, fauna, sentient life) from a world's stored values, with no new column and no epoch: "Microbes: abundant, stable, secure." **Abundance** is the kind's stored 0 to 9 score in the ecologists' DAFOR words (Absent 0, Rare 1 to 2, Occasional 3 to 4, Frequent 5 to 6, Abundant 7 to 8, Dominant 9). **Trend** is Growing, Stable or Declining: the living class is scored again with the star half a billion years on (`main_sequence` brightens, the surface follows the fourth root of the luminosity, the sea is lost past the runaway-greenhouse limit), and a gain or loss of a tenth of its likelihood is a trend; a kind whose milestone (microbes: abiogenesis; flora and fauna: multicellularity; sentient: the technological civilization) came within 0.3 Gyr of the world's age is still Growing. **Security** is how long the class's conditions last, searched in steps of 0.1 Gyr up to 2 Gyr: Dying out (the class scores nothing now, or within 0.1 Gyr), Precarious (within 0.5 Gyr), Marginal (within 2 Gyr, or a score of 2 or less now), Secure. A world with no living class reads Absent for every kind with no trend and no security. The kinds of one world share its living class and so its security; several living classes (GEN.1209) will give one reading per class. The numbers are stand-ins, rows LR1 to LR5 of `stand-in-values.md`. The page and the API use it (UX.251): `GET /api/planets/<id>/life-readings` (and `moons`) works the words out on request from the stored body and its star and stores nothing, and the page shows them in the Life paragraph of At a glance (the kinds that are there, then the classes that suit, no scores) and in a table in the Life section (the words with the score out of 9, and the time the conditions last when it is known).

## The life help pages (DOC.84)

`/help/life` (Help menu, "Life") and the pages under it explain every word the Life section uses, in the pattern of `/help/habitability`: `web/lib/life_help.py` writes the text and `templates/life_help.html` draws it. They are `/help/life` (the four kinds of life and what each score needs), `/help/life/abundance`, `/help/life/trend`, `/help/life/security`, `/help/life/suitability`, `/help/life/stages` (the seven rungs, the three animal grades, the pace of evolution), `/help/life/classes` (the nine classes in one table, with the source tier and adjustment level of each class's air chemistry) and `/help/life/classes/<number>` (one class family: each variant's drinks, eats, breathes and gives off, the liquid, energy, temperature, pressure and supply it needs, how far it can go, its source tier, the words of its stages, and a "Speculative chemistry" note when its air chemistry is at level L3 or L4). Every number on them is read from the code that uses it (`life_readings`, `life_classes`, `life_ledger.json`, `tuning`), and `tests/test_life_help.py` holds each sentence to the style guide and each number to the code. The planet page links to them: each class name goes to `/help/life/classes/<number>#class-<id>`, the Kind, How much, Trend and Security headings carry a "?" link, and the Life section links the stages, suitability and classes. The DOC.57 help framework (a "?" in every page header, a help search, the route-to-help coverage test of DOC.58) is not built yet; these pages are in the Help menu until it is.
