# Where each star type and phenomenon belongs: what is built and what is left (GEN.103)

Boss (2026-10-07 11:47Z) asked for the tiered structure of where stars of different types sit in the galaxy, and for the
kinds of phenomena with similar restrictions, added to the density and type draws. This note is the review of
2026-10-16 turned into a table: which placement rules the code already has, and which three population-dependent rules
for the type draws are still missing, each with a way to calculate it and, for two of them, the constants Boss supplied on 2026-10-16.

Nothing here is a fixed catalogue. Every number a rule uses comes from the sector profile (`galaxy/profile.py`), which is
a function of the galaxy shape, the cell's address and the galaxy seed. A rule below either reads a profile field or
names the physical law that turns profile fields into the number. Where a law needs a constant from the literature the
note says so, marks it for Boss's research, and the code reads it from `tuning.py` with its source named.

Informs: GEN.103, GEN.269, GEN.270, GEN.271, GEN.272, GEN.273, GEN.275.
Status: design, 2026-10-16. Evidence tags: [S] read from the repository, [R] recalled and unconfirmed, [C] computed here.

## 1. What decides where a thing goes

A cell's profile carries nine populations (`profile.populations`, each with a density, a share of the cell's stars, an age
window, a metallicity peak, a velocity dispersion and a rotation rule) and a `kind_factors` map. Two kinds of rule sit
on top of them, and the first is complete.

| Layer | Question it answers | Where it lives | Status |
|---|---|---|---|
| Density | How many systems does the cell hold? | `galaxy/density.relative_density`; `profile.expected_systems` | built |
| Population mix | Which populations are they from? | `galaxy/density.population_densities`; `profile.mix`, `profile.populations` | built (GEN.234) |
| Age and metallicity | How old, how metal rich? | `tuning.STELLAR_POPULATION_AGE_RANGES_GY`; `profile._fe_h`; `bright_stars.py` reads `populations[...].fe_h` | built (GEN.158, GEN.234) |
| Per-kind rate factors | How much of a phenomenon does the place scale? | `profile.kind_factors_at`; `run_sector.py` multiplies `phenomenon_rate_per_star` | built (GEN.236) |
| Scatter placement | Where do the single objects go? | `phenomenon_scatter.py`, one weight per band and longitude bin | built (GEN.278) |
| Type draws by population | Does the place change how many stars a system holds and how far apart they sit? | `generation/multistar.py` | **not built** (section 3) |

## 2. Built: the placement table

Each row is a kind the registry (`generation/phenomena/registry.py`) already places, with the rule that decides where. The
"Follows" column is what the draw is weighted by; no row has a position or a count written into it.

| Kind (registry gate) | Follows | Why | Where it is in the code |
|---|---|---|---|
| star-system (stellar density) | total density; the population draws the age and [Fe/H] | the field stars are the populations | `density.population_densities`, `star_population.py` |
| bright-star (stellar density) | density, population by population | the bright tail of each population's mass function | `bright_stars.py` (pre-placed, per mass cell) |
| black-hole, neutron-star (stellar density) | stars, factor 1 | born where massive stars were, then kicked into the disk's thickness | `kind_factors` 1; `phenomenon_scatter` band density |
| neutron-star, magnetar (young stars) | the young population and its arm contrast, none in the halo or a cluster | too young to have drifted | `SCATTER_CLASS_GATES`, `phenomenon_scatter.py` |
| nebula (hot stars) | the cold gas fraction around the cell's own hot stars | emission and reflection nebulae need a hot star and gas | `kind_factors["nebula"]` |
| molecular-cloud (nebula field) | cold gas fraction times the young share against the solar circle's | the places that form stars | `kind_factors["molecular-cloud"]` |
| supernova-remnant (stellar density) | the gas fraction, none in the halo or in clusters, kept inside a superbubble | the shell shines only where it meets gas | `kind_factors["supernova-remnant"]` |
| planetary-nebula, brown-dwarf (stellar density) | stars, factor 1 | every population makes both | `kind_factors` 1 |
| rogue-planet, comet (stellar density) | the planet occurrence at the cell's [Fe/H], averaged over its populations | planets favour metal-rich stars | `metallicity.planet_occurrence_factor` |
| hypervelocity-star, runaway-star | the central black hole; the stars' own speeds | not a population effect | `wanderer` pass; `flag_fast_stars` |
| quasar | the nucleus | one per galaxy | `nucleus` pass |

The cluster and satellite populations (GEN.160, GEN.237) are placed by the component field, not by a gate: their density
is the components' systems in the cell's cube.

## 3. Left: three rules for the type draws

`multistar-and-compact-systems.md` section 6.6 names three effects of the population on a system's stars. None is in
`multistar.py` yet: the star count, the pairing and the axis limit read the star's mass and nothing about the place.
Each rule below is a function of fields the profile already holds. The constants are the only things to look up, and
they sit in `tuning.py` with a source.

### 3.1 Close pairs rise as metallicity falls

Moe, Kratter and Badenes (2019, ApJ 875, 61; arXiv 1808.02116) corrected five surveys for their selection biases and found
one consistent trend: the intrinsic fraction of solar-type primaries (about 0.6 to 1.5 solar masses) with a companion inside
10 AU (periods under 10^4 days) is **53% +-12% at [Fe/H] -3.0, 40% +-6% at -1.0, 24% +-4% at -0.2 and 10% +-3% at +0.5**
[S: the abstract, read through search results; the paper's full text could not be opened from this environment].
Badenes et al. (2018, ApJ 854, 147) found the same direction in APOGEE: metal-poor stars ([Fe/H] below about -0.5) are
two to three times as often in multiples as metal-rich ones (above about 0). Bate (2019) reproduces it in
radiation-hydrodynamic simulations, where low opacity lets the disc fragment more and the companions migrate inside
10 AU [R].

Rule: the chance that a drawn pair is a close one (an axis under 10 AU) scales with the star's [Fe/H], which the cell
already draws per population (`PopulationProfile.fe_h`, peak and scatter). The existing multiplicity tables
(`BINARY_SYSTEM_PROBABILITY_BY_SPECTRAL_CLASS`, `STAR_COUNT_PROBABILITIES`) keep their mass dependence; the metallicity
enters as a factor `F(fe_h) / F(0)`, where `F` is the straight line between the four points above (in [Fe/H]), held flat
past the ends. `F(0)` is about 20%, so a star at -1.0 gets about twice the close pairs of a solar one and a star at +0.5 about
half. The factor is a function of the star's own [Fe/H] and nothing else, so no population is named.

What the sources limit:

- **Which stars the four points describe.** Field solar-type dwarfs and subgiants, pooled over the disc and halo samples the
  five surveys hold; they are not four populations. The -3.0 end rests on the Carney-Latham halo sample and the +0.5 end on
  Kepler and APOGEE, so the line between them is an interpolation of a pooled trend, not a measurement per population.
- **Mass.** The points are for 0.6 to 1.5 solar masses. Moe et al. argue the close fraction of massive stars does not change
  with metallicity, so the factor should fade to 1 for primaries above about 2 solar masses [R]. Boss's text gave
  0.6 to 1.2 solar masses and "periods below 10^3.5 days"; 10 AU around a solar-mass star is a period of about 10^4.06 days, so
  the paper's 10^4 days is the consistent form.
- **Wide pairs.** Boss's text says the dependence vanishes beyond about 200 AU. El-Badry and Rix (2019, MNRAS Letters; arXiv 1809.06860)
  say the wide fraction is flat in [Fe/H] beyond about 250 AU and anticorrelated below about 200 AU, with the fraction at [Fe/H] -1 about three times
  that at +0.5 for 50 to 100 AU. Two later Gaia studies disagree with "flat": Hwang et al. 2021 (arXiv 2010.02920) find the wide fraction peaks
  near [Fe/H] 0 and falls to both sides for 10^3 to 10^4 AU, and Niu et al. 2022 (arXiv 2205.04444) say the difference comes from the
  separation range and the stellar temperature, with the sign flipping between 4000 and 7500 K. The rule therefore moves the
  close fraction only and leaves the wide fraction alone (the choice Boss's text states) until he decides.
- **Iron or alpha.** Mazzola et al. (2020, MNRAS 499, 1607; arXiv 2007.09059) find the close fraction follows the alpha abundance more
  steeply than iron, with a high alpha abundance suppressing it at most [Fe/H]. The generator has no alpha abundance. Used with [Fe/H] alone the rule would give the thick disc
  more close pairs than a high-alpha sample shows. [R: from the abstract only]

### 3.2 Massive multiples sit in young regions

Observed: high-mass stars are far more often in multiples than low-mass ones, and the most massive sit in clusters and
associations. The mass dependence is already in `BINARY_SYSTEM_PROBABILITY_BY_SPECTRAL_CLASS` and
`STAR_COUNT_PROBABILITIES`; the place dependence is not. The `prevalence` argument of `draw_member_tree` already scales
the chance of five or more stars, so the place only has to set it.

Rule: `prevalence` is the young share against the solar circle's, clamped to 0 to some cap, where the young share is the profile's
`populations["young"].share` over the reference `profile._reference_young_share(shape)` that `kind_factors_at` already
uses. That makes the knob 1 where a cell looks like the solar circle, above 1 in arms and associations and 0 in the halo. The cap on the rise is the one constant (a design default for Boss, since the five-or-more column is a tail).
The five-or-more columns are the ones that scale, which is how the knob works today. 
### 3.3 Wide pairs survive less where it is dense

Two limits cut the widest axis a pair can keep. Both are functions of the cell and the second also of the stars' age.

**The encounter limit.** Passing stars and clouds unbind soft pairs (Heggie 1975; Bahcall, Hut and Tremaine 1985). Jiang and
Tremaine (2010, MNRAS 401, 977; arXiv 0907.2952) restate the 1985 result: the semi-major axis at which **half of
the pairs of age t survive** is

    a_1/2 = 0.002 * (M_1 + M_2) * sigma / (G * rho_2 * t)

with sigma the velocity dispersion of the field, t the pair's age and rho_2 the density of the perturbers weighted by
their mass (units of solar masses squared per cubic parsec, about the mass density times the mean perturber mass)
[S: a secondary source quoting the paper; the constant 0.002 is that source's, and what it folds in (the Coulomb logarithm, the encounter integral) was not checked].
The same source gives its worked value: 3.1 x 10^4 AU for M_1 + M_2 = 2, sigma = 50 km/s, rho_2 = 0.03 and t = 10 Gyr, which I
reproduced to 31,300 AU [C]. It is linear in 1/rho_2 and 1/t, which is what the code needs: a_1/2 follows the cell's density
and the star's age with nothing stored by place. In the solar neighbourhood (sigma 30 km/s, rho 0.08 solar masses per pc^3, a perturber
mass of 0.4) it gives 1.8 x 10^4 AU, near Boss's 2 to 5 x 10^4 AU and Jiang and Tremaine's 50% survival near 30,000 AU.

**The tidal limit.** A pair of total mass M has a Jacobi radius of (G M / (4 Omega^2 - kappa^2))^(1/3), and
4 Omega^2 - kappa^2 = 4 A (A - B) from Oort's constants. Boss's text wrote the last form with 2 A (A - B), which does not match
its own first form; with A = 15.3 and B = -11.9 km/s/kpc the correct one gives 4 Omega^2 - kappa^2 = 1.75 x 10^-30 s^-2 (his
"about 2 x 10^-30") and **1.73 pc** for 2 solar masses, which matches Jiang and Tremaine's 1.7 pc; the 2 A (A - B) form would give 2.18 pc [C].
Section 6.6 of [multistar-and-compact-systems.md](multistar-and-compact-systems.md) takes this form from the galactic potential at the cell,
the one the orbital update uses. Where the cell holds only its density, (M / (3 rho))^(1/3) is the same limit written with the local
mass density (1.9 pc for 0.1 solar masses per pc^3 [C]) and serves as the check.

The widest axis a wide pair may keep is the smaller of the two. It takes the place of the fixed 10,000 AU cap
(`WIDE_BINARY_SEPARATION_MAX_AU`), as section 6.6 now says. Today's shape, a 2 solar mass pair of age 10 Gyr, the density of the
model's `relative_density` times `E / edge^3` times `tuning.MEAN_SYSTEM_MASS_SOLAR` (0.4), perturber mass 0.4 [C]:

| Place | rho (solar masses per pc^3) | sigma (km/s) | Encounter limit a_1/2 | Tidal limit | Widest axis kept |
|---|---|---|---|---|---|
| 8.2 kpc out, in the plane | 0.048 | 30 | 29,000 AU | 500,000 AU | 29,000 AU |
| bulge, 500 pc from the centre | 6.9 | 120 | 800 AU | 95,000 AU | 800 AU |
| bulge, 100 pc from the centre | 77 | 120 | 73 AU | 42,000 AU | 73 AU |
| 10 pc from the centre | 660 | 120 | 9 AU | 21,000 AU | none (under the 50 AU lower edge) |
| 2 pc from the centre | 61,000 | 150 | under 1 AU | 4,600 AU | none |

So the encounter limit, not the tide, decides everywhere but the quiet disc, and the old fixed cap of 10,000 AU is exceeded
in the disc (29,000 AU) and far too wide in the bulge. That is the intended consequence of 6.6, and it is a visible one: a
bulge sector would hold almost no wide pairs. The reference density `rho_2` and sigma come from the profile (`gas.mass_solar`
and `dark_matter` are left out of rho, which would shorten the limit a little). The sigma used is the share-weighted rms of the
cell's populations (`tuning.STELLAR_POPULATION_DISPERSION_KMS`). The age is the star's own, so a young pair keeps a wider axis than an old one.

## 4. How the rules reach the draw

`run_sector.py` already passes `profile.kind_factors` into the fill. The three rules need the profile's populations and
density in the system draw, so the system config gets one more input, a small `PlaceTerms` record built once per cell
from the profile (the [Fe/H] per population, the young-share prevalence, the density in solar masses per cubic parsec),
read by `multistar.draw_member_tree` (section 3.2 through `prevalence`), the pair-pick (3.1) and `wide_cap` (3.3). It is
a function of the cell, so the same seed gives the same stars wherever a run lands (the GEN.39 rule). It stays out of
`tuning.py`'s stellar evolution constants, so Core 3's evolution work is not touched.

## 5. What closes GEN.103

- 3.1 and 3.3 have their constants (the four points; the 0.002 coefficient). 3.2 needs a design cap on the prevalence.
- The choices Boss still has to make are in section 6.
- GEN.269 and GEN.270 to GEN.275 are the multi-star items that own the code (Core 3), so this note is their input: the
  build goes there, and GEN.103 closes when they have taken the three rules in, or the rules are written as items of
  their own.

## 6. Check of the pasted reference list (2026-10-16, 6:15 PM EDT)

Boss pasted the research; each entry was checked against what search results show of the source. The papers themselves
could not be opened (arxiv.org is not reachable from this environment), so "confirmed" means the title, authors, journal and the quoted
result appear together in at least one source result; numbers that rest on a secondary source are marked.

| Pasted citation | Result | What was found |
|---|---|---|
| Moe, Kratter and Badenes 2019, ApJ 875, 61 | confirmed | Title "The Close Binary Fraction of Solar-type Stars Is Strongly Anticorrelated with Metallicity", arXiv 1808.02116. Points 53, 40, 24, 10% (section 3.1). Periods under 10^4 d and 0.6 to 1.5 solar masses, not 10^3.5 d and 0.6 to 1.2. |
| Badenes et al. 2018, ApJ 854, 147 | corrected | Volume and page right. The title is "Stellar Multiplicity Meets Stellar Evolution and Metallicity: The APOGEE View", with Mazzola among the authors; the pasted title, "Stellar multiplicity in the APOGEE survey. I. ...", and the author list (Maoz, Ciotti, Krajnovic, Shen, Simon) are wrong. Finding: metal-poor stars have 2 to 3 times the multiplicity of metal-rich. |
| El-Badry and Rix 2018, MNRAS 480, 4884 | corrected | The metallicity paper is El-Badry and Rix 2019, MNRAS Letters, "The wide binary fraction of solar-type stars: emergence of metallicity dependence at a < 200 AU" (arXiv 1809.06860). MNRAS 480, 4884 (2018) (DOI sty2186) is a different 2018 Gaia wide-binary paper; its authors were not checked. The pasted title is wrong. |
| El-Badry, Rix and Heintz 2021, MNRAS 506, 2269 | corrected | Journal right. Title is "A million binaries from Gaia eDR3: sample selection and validation of Gaia parallax uncertainties" (arXiv 2101.05282), not "...validation of theoretical models". Its sample runs to 1 pc. A truncation of the separation distribution at 2 to 5 x 10^4 AU was not found in it: one review says chance alignments dominate beyond about 30,000 AU, and a different study finds a break near 0.1 pc. Could not confirm. |
| Kratter and Lodato 2016, ARA&A 54, 271 | confirmed | "Gravitational Instabilities in Circumstellar Disks", pages 271-311, arXiv 1603.01280. |
| Bate 2018, MNRAS 484, 2341 | corrected | The metallicity work is Bate (2019), "The statistical properties of stars and their dependence on metallicity" (arXiv 1901.03713; the 2014 "effects of opacity" paper is MNRAS 442, 285). It found low metallicity gives more small-scale fragmentation, matching the close-binary trend. The pasted year and title are wrong, and 484, 2341 could not be confirmed. |
| Heggie 1975, MNRAS 173, 729 | confirmed in substance | The result ("soft binaries widen and dissolve, hard ones harden") is cited as Heggie's law in the search results; the volume and page were not seen. |
| Bahcall, Hut and Tremaine 1985, ApJ 290, 15 | corrected | Volume and page right. The title is "Maximum mass of objects that constitute unseen disk material", not "Maximum sizes of compact perturbations in the Galaxy". |
| Jiang and Tremaine 2010, MNRAS 401, 977 | confirmed | "The Evolution of Wide Binary Stars", arXiv 0907.2952. 1.7 pc Jacobi radius for a solar-mass pair confirmed; 50% survival near 30,000 AU is a secondary source's reading. |
| Stegmann et al. 2024 | corrected | Astrophysical Journal Letters 972, L19, "Close Encounters of Wide Binaries Induced by the Galactic Tide: Implications for Stellar Mergers and Gravitational-wave Sources" (arXiv 2405.02912), by Stegmann, Vigna-Gomez, Rantala, Wagg, Zwick, Renzo, van Son, de Mink and White; not Physical Review D, and Antonini and Perets are not authors. It finds pairs wider than about 10^5 AU are torn apart by the tide, which does not support "natal cluster dissolution". |
| Gammie criterion, t_cool under 3 / Omega | not checked | Standard disc-fragmentation criterion; no source needed for the rule. |
| Boss's diffusive formula | wrong as typed | a = (G M^2 / (8 pi G^2 n M_*^2 ln Lambda T))^(1/2) is not a length (units of mass times time) and holds no velocity dispersion; replaced by the form in 3.3. |
| Jacobi form "2 A (A - B)" | wrong as typed | The first form, 4 Omega^2 - kappa^2, equals 4 A (A - B) (section 3.3). |

## 7. Choices for Boss

1. **Which population the four close-binary points describe.** (a) Use them as the pooled field trend, a function of the star's own
   [Fe/H] (recommended: it is what the paper measures, and needs no population). (b) Split by population with different lines: no
   source gives that. (c) Add an alpha-abundance term (Mazzola et al. 2020): needs an alpha abundance the generator does not draw.
2. **Mass range of the metallicity rule.** (a) Fade the factor to 1 above about 2 solar masses (recommended: the paper says massive
   stars do not change). (b) Apply it to every mass.
3. **The wide-pair metallicity dependence.** (a) None (recommended: it is Boss's stated rule and El-Badry and Rix agree above 250 AU). (b) The
   Hwang et al. non-monotonic curve for 10^3 to 10^4 AU: a table for a different separation range, not yet read.
4. **The form of the encounter limit.** (a) A sharp cut at a_1/2 (recommended: simplest, deterministic, one coefficient).
   (b) A survival probability falling through a_1/2 (the papers say the real curve is gradual): needs the survival curve from
   Jiang and Tremaine's figures, which could not be read.
5. **The fixed 10,000 AU cap.** The derived limit replaces it (6.6), so the quiet disc gets wider pairs (29,000 AU) and a bulge
   sector almost none. (a) Accept (recommended: it is what 6.6 says). (b) Keep 10,000 AU as an upper bound as well.
