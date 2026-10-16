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
Status: 3.1 and 3.3 built (GEN.103); 3.2 left; check TEST.196 queued. Design, 2026-10-16. Evidence tags: [S] read from the repository, [R] recalled and unconfirmed, [C] computed here.

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

`multistar-and-compact-systems.md` section 6.6 names three effects of the population on a system's stars. 3.1 and 3.3
are built (GEN.103: `generation/place_terms.py`, `StarSystem._companion`, `wide_binary.py`); 3.2 is not, and the scatter's
member trees (`multistar.py`) read the star's mass and nothing about the place.
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

Rule (Boss, 6:12 to 6:15 PM EDT, decided): the close-companion chance follows one pooled line in [Fe/H], read at the star's own
[Fe/H] (`PopulationProfile.fe_h`, peak and scatter, drawn per population by the cell), for **every primary mass, with no fade**.
`F` is the straight line between the four points above, held flat past the ends; `F(0)` is 20%. The factor
`F(fe_h) / F(0)` multiplies the close share the generator already has and nothing else moves:

- Today a pair is drawn in two rolls: `BINARY_SYSTEM_PROBABILITY_BY_SPECTRAL_CLASS[class]` for "has a companion" (O 0.90, B 0.65, A 0.55,
  F 0.47, G 0.44, K 0.40, M 0.26), then `WIDE_BINARY_DEFAULT_CHANCE` (0.5) for wide against close (`system.py`, `_should_generate_binary`
  and the split after it). So today's close share is `p_close0 = P_binary(class) * (1 - 0.5)` (G 0.22, O 0.45) and the wide share
  `p_wide = P_binary(class) * 0.5`. The generator's close pairs sit at 0.05 to 0.25 AU (`binary.py`), all inside the 10 AU the points describe.
- New: `p_close = min(1, p_close0 * F(fe_h) / F(0))`; `p_wide` is unchanged (answer 3: no dependence); one roll then picks close, wide or
  single from `p_close`, `p_wide` and the rest, with `p_close + p_wide` held to at most 1: where the sum would pass 1 the wide share is cut to `1 - p_close`. That happens
  for O stars below [Fe/H] about -0.25 and B stars below about -1.1 (the pair table is high enough that the companion chance saturates), and it is the one place the wide share
  moves with [Fe/H], an exception to answer 3 that the wide-flatness check skips.
- A G dwarf (p_close0 = 0.22) then gets 0.58, 0.44, 0.26 and 0.11 at [Fe/H] -3.0, -1.0, -0.2 and +0.5, against 53, 40, 24 and 10%
  (+-12, 6, 4, 3) in the paper [C]: the generator's own table, with the halved wide split, already reproduces the solar-neighbourhood
  value, so no recalibration is needed for the solar-type bin.

The factor is a function of the star's own [Fe/H] and nothing else, so no population is named.

What the sources limit:

- **Which stars the four points describe.** Field solar-type dwarfs and subgiants, pooled over the disc and halo samples the
  five surveys hold; they are not four populations. The -3.0 end rests on the Carney-Latham halo sample and the +0.5 end on
  Kepler and APOGEE, so the line between them is an interpolation of a pooled trend, not a measurement per population. Boss chose
  the pooled line (answer 1) over a per-population one, which no source gives.
- **Mass: all masses, which the sources do not support above about 1.5 solar masses.** The points are for 0.6 to 1.5 solar masses.
  Moe, Kratter and Badenes argue the close fraction of massive stars does not change with metallicity [R]; for O stars and early-B stars
  Moe and Di Stefano (2013, eclipsing fractions in the Small and Large Magellanic Clouds and the Milky Way), Dunstall et al. (2015) and
  Almeida et al. (2017) find the close fraction in the Large Magellanic Cloud ([Fe/H] about -0.4) consistent with the Milky Way's
  [S: search summaries; the papers could not be opened]; the BLOeM survey (arXiv 2503.21936) finds early-B stars at the Small Cloud's
  metallicity ([Fe/H] about -0.7) *more* often close than in the Milky Way, so early B may follow the trend while O stars do not.
  Boss's decision is to apply the line to every mass and let the check (3.1b) show where the numbers leave the published ranges. What the
  line does to a massive star is in the table below: it **will** push O stars below [Fe/H] about -0.7 past the observed
  ceiling, so the check is expected to flag those bins on its first run.
- **Wide pairs.** Answer 3: no metallicity dependence, so the rule moves the close share only. El-Badry and Rix (2019, MNRAS Letters; arXiv 1809.06860)
  say the wide fraction is flat in [Fe/H] beyond about 250 AU and anticorrelated below about 200 AU; Hwang et al. 2021 (arXiv 2010.02920) find it peaks
  near [Fe/H] 0 and falls to both sides for 10^3 to 10^4 AU, and Niu et al. 2022 (arXiv 2205.04444) say the sign depends on separation range
  and temperature. The generator's wide pairs start at 50 AU (`WIDE_BINARY_SEPARATION_MIN_AU`), so the choice leaves the 50 to 200 AU
  anticorrelation out; that is Boss's rule.
- **Iron or alpha.** Mazzola et al. (2020, MNRAS 499, 1607; arXiv 2007.09059) find the close fraction follows the alpha abundance more
  steeply than iron, with a high alpha abundance suppressing it at most [Fe/H]. The generator has no alpha abundance. Used with [Fe/H] alone the rule would give the thick disc
  more close pairs than a high-alpha sample shows. [R: from the abstract only]

What the pooled line does to each class (close share `min(1, p_close0 * F/F(0))`, from the generator's table, all [C]):

| [Fe/H] | F | F / F(0) | O | B | A | F | G | K | M |
|---|---|---|---|---|---|---|---|---|---|
| -3.0 | 53% | 2.65 | 1.00 | 0.86 | 0.73 | 0.62 | 0.58 | 0.53 | 0.34 |
| -1.0 | 40% | 2.00 | 0.90 | 0.65 | 0.55 | 0.47 | 0.44 | 0.40 | 0.26 |
| -0.7 | 34% | 1.70 | 0.77 | 0.55 | 0.47 | 0.40 | 0.37 | 0.34 | 0.22 |
| -0.4 | 28% | 1.40 | 0.63 | 0.45 | 0.39 | 0.33 | 0.31 | 0.28 | 0.18 |
| -0.2 | 24% | 1.20 | 0.54 | 0.39 | 0.33 | 0.28 | 0.26 | 0.24 | 0.16 |
| 0.0 | 20% | 1.00 | 0.45 | 0.33 | 0.28 | 0.23 | 0.22 | 0.20 | 0.13 |
| +0.5 | 10% | 0.50 | 0.23 | 0.16 | 0.14 | 0.12 | 0.11 | 0.10 | 0.07 |

### 3.1b The check that the numbers act as expected (TEST.196)

Boss (6:15 PM EDT): "build into that a check to see if the numbers act as expected." It is TEST.196, owner Core 4, in the census harness next
to TEST.186 (Core 3's; ask before editing `census.py`). It fills seeded sectors across the [Fe/H] range, counts the share of
primaries with a companion inside 10 AU by primary-mass bin and [Fe/H] bin, and compares each bin to the expected range below plus
the binomial allowance of its sample (three standard errors), and shows the wide share is flat in [Fe/H]. It fails the run when a bin with a published
range is outside it, prints every bin (also the ones with no range), and writes the table to the audit report. The ranges are reference values for the check
and are never read by the draw.

| Primary mass bin | Expected close share (inside 10 AU) | Source | What the check does |
|---|---|---|---|
| 0.6 to 1.5 Msun (G, late F, early K) | 53 +-12% at [Fe/H] -3.0, 40 +-6% at -1.0, 24 +-4% at -0.2, 10 +-3% at +0.5 | Moe, Kratter and Badenes 2019 [S: abstract] | fails when outside the range widened by the sampling error |
| 7 +-2 Msun (mid B) | at least 23 +-7% (companions at periods 100 to 1500 d only, mass ratio 0.1 to 1); [Fe/H] about 0 | Moe and Di Stefano 2017 [S: search summary of the abstract] | fails below the lower bound at [Fe/H] -0.2 to 0 |
| 28 +-8 Msun (O) | between 35 +-12% (periods 100 to 1500 d) and about 69 +-9% (Sana et al. 2012, periods to about 3000 d, intrinsic binary fraction) | Moe and Di Stefano 2017; Sana et al. 2012 [S: search summaries] | fails outside the span at [Fe/H] -0.2 to 0 |
| O and early B across [Fe/H] -0.7 to +0.1 | the same value within about 30% (no trend) | Moe and Di Stefano 2013; Dunstall 2015; Almeida 2017 [S] | fails when the ratio of the lowest to the highest [Fe/H] bin leaves 0.7 to 1.3 (**expected to fail**, see below) |
| below 0.6 Msun, 1.5 to 5 Msun | no number verified in this environment; the close fraction rises with mass | Moe and Di Stefano 2017 (not opened) | printed, never fails, until Boss does the research |
| wide share (50 AU and out), all masses | no change with [Fe/H] (answer 3) | the rule itself | fails when two [Fe/H] bins differ by more than the allowance; skips the O and B bins where the companion chance saturates (see above) |

Predicted first-run result, from the table above: the G-dwarf bins pass at all four points; the O bin at [Fe/H] -0.7 and below comes out at 0.77 to 1.00,
above the Sana ceiling, and the O and early-B flatness test fails (factor 1.7 between -0.7 and 0); at +0.5 the O share falls to 0.23, on the lower edge of the 35 +-12% reference
(35 - 12 = 23). So the safeguard is built to fire on the massive stars: it is the evidence for the next decision, not a surprise.
If it fires, the options are (a) cap the factor at the observed 30% for primaries above about 2 solar masses (the fade Boss declined), (b) accept it and loosen the
check, (c) fade only the O class. That is Boss's call when he sees the report. Until he picks, TEST.196 lands with the massive-star flatness test marked as an expected
failure, so the run shows it every time without hiding the other bins.

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

The widest axis a wide pair may keep is the smaller of the two (answer 4: a **sharp cut** at a_1/2, a pair wider than that does not exist). It takes the
place of the fixed 10,000 AU cap (`WIDE_BINARY_SEPARATION_MAX_AU`, answer 5), as section 6.6 now says: the axis is still drawn log-uniform
from the 50 AU lower edge, up to the limit instead of 10,000 AU. A cell whose limit is under 50 AU keeps no wide pair: the wide roll then gives the primary alone
(the pair dissolved), which is the whole effect in the bulge. The wide *chance* (0.5 of the binary chance) does not change where the limit is
above 50 AU; scaling it by the share of the log range kept would need a natal upper bound the sources do not give, so it is left out (a design default). Today's shape, a 2 solar mass pair of age 10 Gyr, the density of the
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

Built. The close-pair factor reads the primary's own [Fe/H] (`primary_star.metallicity_fe_h`), so it needs nothing from the cell.
The wide-pair limit needs the cell's density and dispersion: `place_terms.from_profile(profile)` builds a `PlaceTerms` once per
cell (`FillContext.place_terms`), `FillContext.apply` and the pre-placed bright stars put it on `SystemConfig.PLACE_TERMS`, and a
system with none takes the Sun's neighbourhood (`tuning.SOLAR_NEIGHBOURHOOD_*`). `StarSystem._companion` rolls none, close or wide
in one draw and `WideBinaryPair` draws its axis up to `wide_pair_limit_au`. It is a function of the cell, so the same seed gives
the same stars wherever a run lands (the GEN.39 rule), and it stays out of the stellar evolution constants. Section 3.2 would need
the same record to carry the young-share prevalence.

## 5. What closes GEN.103

- 3.1 and 3.3 are built to Boss's five answers (section 7). 3.2 needs a design cap on the prevalence and belongs with the multi-star
  items (GEN.269 and GEN.270 to GEN.275, Core 3), which can read `PlaceTerms` once they need it.
- TEST.196 (3.1b) follows, in the census harness (TEST.186, Core 3's).

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

Sources added for the TEST.196 ranges (2026-10-16, 6:28 PM EDT; same method: search summaries only, the papers could not be opened):

| Citation | Result | What was found |
|---|---|---|
| Moe and Di Stefano 2017, ApJS 230, 15 | confirmed | "Mind Your Ps and Qs" (arXiv 1606.05347). Binary fraction for periods 100 to 1500 d and mass ratio 0.1 to 1: 23 +-7% for mid-B primaries (7 +-2 Msun), 35 +-12% for O primaries (28 +-8 Msun). The companion count rises from 0.48 +-0.04 (solar-type) to 2.1 +-0.3 (O), but that counts all separations. Table 13 holds the per-mass multiplicity statistics; its close-fraction rows could not be read. |
| Sana et al. 2012, Science 337, 444 | confirmed | Intrinsic binary fraction 0.69 +-0.09 for Galactic O stars, periods up to about 3000 d. |
| Sana et al. 2013 (VLT-FLAMES Tarantula Survey VIII, A&A) | confirmed in substance | 0.51 +-0.04 for log P 0.15 to 3.5 in 30 Doradus; volume and page not seen. |
| Moe and Di Stefano 2013 | confirmed in substance | Eclipsing fractions of early-B stars 0.70%, 0.69% and 1.00% (Small Cloud, Large Cloud, Milky Way), read as no metallicity dependence; title and journal not seen. |
| Dunstall et al. 2015; Almeida et al. 2017 | confirmed in substance | Large Cloud early-B and O close fractions consistent with the Milky Way's after completeness correction; journals not seen. |
| BLOeM (arXiv 2503.21936) | confirmed in substance | Early-B close fraction at the Small Cloud's metallicity higher than the Milky Way's and the Large Cloud's, hinting at an anticorrelation for B stars; O stars constant (observed 45 +-4%, modelled about 0.70). |
| Kobulnicky et al. 2014 ("55% for periods under 5000 d") | not confirmed | The attribution did not survive the search; not used. |
| Moe and Di Stefano 2013, "frequency of very close companions scales as M^0.7" | quoted by a secondary source | Not used as a range. |

## 7. Boss's answers (6:12 to 6:15 PM EDT, 2026-10-16)

1. **Close-binary trend.** One pooled field line (53, 40, 24, 10%), read at the star's own [Fe/H] (3.1).
2. **Mass range.** Every star mass, no fade, with a check that the numbers act as expected: TEST.196 (3.1b). The note records that the sources
   support the line only for 0.6 to 1.5 solar masses, and that the check is expected to flag the massive bins.
3. **Wide pairs.** No metallicity dependence.
4. **Encounter limit.** A sharp cut at a_1/2.
5. **The 10,000 AU cap.** Replaced by the derived limit (3.3); a cell whose limit is under 50 AU keeps no wide pair.

Design defaults chosen here, for Boss to overrule: the companion roll saturates at 1 and cuts the wide share (3.1); the wide chance is not scaled by the
log range kept (3.3); a dissolved wide pair leaves the primary alone (3.3).
