# The life timeline and the oxygen life makes (GEN.750, GEN.855)

Status: built with GEN.855 and GEN.750 (one change, one epoch). Code: `physics/life_timeline.py`, `physics/biogenic_oxygen.py`, `generation/life.py` (`_draw_schedule`, `_build_oxygen`), `generation/evolution.py`. Constants: `tuning.LIFE_*` and `tuning.SUNLIT_OCEAN`. Every value that is a stand-in is in the registry (`docs/design/research/gap-filling-plan.md`, section "Stand-ins of the life timeline and the biogenic oxygen").

## 1. What was wrong

Before this change a world's stage of life came from a table of paces per chemistry, and oxygen came from an atmosphere draw that knew nothing of life. A world could have complex life with no oxygen in its air (the population tests found no technological civilization at all), and a world with oxygenic photosynthesis stayed anoxic. Boss's bug GEN.855: the air of a living world must hold the oxygen its life made, by the time its life has had to make it.

## 2. The timeline (GEN.750)

Five stages, in the order of `tuning.LIFE_MILESTONES`: abiogenesis, photosynthesis (the first harvesting of light), complex cells, multicellularity, technological civilization. Each is a hard step in the sense of Carter (1983) and Watson (2008): a rare event with an exponential waiting time whose mean (`LIFE_STEP_MEANS_GYR` = 0.5, 1.5, 2.5, 3, 60 Gyr, times `LIFE_STEP_SCALE` = 1; Boss chose them, 9:03 PM EDT 10-28) is often longer than the world is old. The wait for a step starts when the one before it is done. A world has reached every stage whose age is less than the age of its star. No Earth history is copied.

Three things shape a world's own clock:

* **Temperature.** A Boltzmann factor on the rate of every step, `exp(E / k (1/T_ref - 1/T))`, with E = 0.30 eV (`LIFE_ACTIVATION_ENERGY_EV`) and T_ref = 288 K, kept inside 0.25 to 2.5. The pace words ("fast", "normal", "slow") name that factor. Metabolic theory gives 0.65 eV for a cell's chemistry; a hard step is a rare mutation and selection event whose dependence the record cannot pin, and 0.65 eV made a warm world three times faster than the reference.
* **The kind of light harvesting.** The second step is scaled by the chemistry's `metabolism` (`LIFE_FIRST_LIGHT_SCALE`): Retinal 0.4, anoxygenic photosynthesis 0.6, radiotrophy 1.0, oxygenic photosynthesis 1.0. It keeps the order of the record: a pigment and a proton gradient first, one reaction centre next, the water-splitting machinery last.
* **Oxygen.** Complex cells need free oxygen, so their wait starts only when the air holds some: at once if it already does (an abiotic oxygen air, 200 Pa or more), else a delay after the first photosynthesis when the life is oxygenic and the world can make the oxygen (section 3), and never otherwise. A world whose life is not oxygenic and whose air has no oxygen therefore stops at photosynthesis, which is the answer the model gives for Retinal, anoxygenic and radiotrophic life.

A requested technological civilization (`INTELLIGENT_LIFE`) is today a timeline conditioned on finishing inside the world's age: `forced_waits` draws timelines on a stream of its own until one fits and squeezes the last one to fit when none does. **Boss ruled at 5:08 PM EDT 2026-10-30 that a civilization is never forced**, so a request must only tweak probabilities; the existing forcing of the `+intelligent_life` option is a question for Boss, and the retry item GEN.913 was withdrawn. The run's stream moves by one draw whatever happens.

Measured shares (40,000 draws, star age uniform 3 to 13 Gyr, oxygenic life, reference temperature, oxygen 2 Gyr after the first photosynthesis): complex cells 65 in 100, animals 41 in 100, a civilization 2 in 100 (for G, K and M stars of 1 to 10 or 13.8 Gyr: animals 21, 33, 23 and a civilization 0.6, 1.5, 0.7). The scale (1.0) is the one knob Boss can retune.

## 3. The oxygen (GEN.855)

Oxygenic life (Chlorophyll a and Blue-Optimized Porphyrins; `metabolism` is `oxygenic`) needs a sunlit sea: `hydrosphere` is a surface ocean or a hycean ocean (`SUNLIT_OCEAN`). The oxygen of the air is the balance of a source and its sinks, in moles of O2 a year:

* **Source S**, the organic carbon its life buries: planet area, the ocean share, the square root of the photosynthetically active light, a Boltzmann factor of 0.32 eV, to zero outside the liquid range of the sea.
* **Reductant flux R**: volcanic and metamorphic gases that take oxygen as they come out; it goes with the heat the interior gives off and the redox of the mantle.
* **Oxidative weathering k N**: rocks take oxygen in proportion to what the air holds.
* **A store B** of reductants (ocean iron, crust sulfides, and the reduced gases of the air itself) that takes the first oxygen, so the air holds none until B / (S - R) has passed (the Great Oxidation came 0.3 to 1 Gyr after oxygenic photosynthesis began on Earth).

Once the store is spent, N settles to (S - R) / k, a time short beside the Gyr of the schedule, within the fire limit (35 percent of the molecules). The numbers are Earth's, used as the test point and the scale, not as a preset: Earth's light, water, heat and redox reproduce 22.7 kPa and a delay of 0.72 Gyr, and another world gets another oxygen. A hycean world burns away its hydrogen first (that is its store), so its oxygen comes late.

The air is bookkept in moles: n_i = x_i P A / (g M), the new pressure is g sum(n_i M_i) / A, and partial pressures follow from the mole shares. The reduced gases (H2, CO, CH4, NH3, H2S, C2H6) burn to their oxides, and the water goes to the ocean. After the oxygen is added everything that reads the air is made again: the greenhouse warming and surface temperature, the stored air, ozone and the radiation doses (on a stream of the body's own), the colour of the atmosphere, and the habitability score. A world that already holds the oxygen is left as it is.

## 4. What the model leaves out

Loss of oxygen to space, the hydrogen escape that adds oxygen, the phosphorus feedback and land plants (they would change S by a factor inside the 1.7 to 4 of this model's tweak level). Oxygenic chemistries are not viable under an ice lid here (the light a lid lets through is left out), so an ice-covered ocean makes no oxygen; whether life under ice could, is an open question.

## 5. What a player sees

A habitable world with oxygenic life older than the delay has an air with oxygen in the range the animals need (8 to 50 kPa for the habitability score `phi_cpx`), and the sheet shows the stage the schedule reached. Worlds with other light harvesting show no oxygen and stop at photosynthesis. The numbers are the same on every run of the same seed.
