# Multi-star systems at scatter time: census, limits, probabilities and the stored tree

Boss (2026-10-15, 5:51 PM EDT): probabilities accurate but "no forced representation in a sector (except if one landed there in a scatter before generation)"; stars and star-like objects (former stars, protostars) are generated, and "if a star, black hole, neutron star, etc has a companion I want it generated at scatter time", including neutron stars with neutron stars, black holes with neutron stars, black holes with black holes and stars with neutron stars or black holes, because "one of the points of our scattering is to represent gravity and such systems would definitely impact gravity". He asked for research into how many multi-star (or star-like object) stable systems exist, the theoretical limits, and the probability of each, then an expansion of our multi-star definitions.

This note does that. It extends [multistar-and-compact-systems.md](multistar-and-compact-systems.md) (the stability criteria, the 2 to 7 star hierarchy, the compact-primary rules; GEN.128 to GEN.130), [phenomena-catalogue-and-fill.md](phenomena-catalogue-and-fill.md) (the kind list and the fill stages; GEN.250, GEN.251) and [exotic-environments-planets-and-compact-binaries.md](exotic-environments-planets-and-compact-binaries.md) (Peters time, V404 Cygni, survival of wide companions). It does not repeat their numbers; it says where this note changes them (section 10).

Informs: GEN.129, GEN.130, GEN.223, GEN.109, GEN.250, GEN.251, GEN.100, GEN.242, GEN.243, GEN.252, and the new items in section 11.

Status: design, 2026-10-15; nothing here is built. Boss approved the plan at 6:08 PM EDT and widened it to every observed and theoretical type; section 13 onward is that catalogue. Where Boss's words leave a gap the note takes a default, marked "Default" and listed in section 12.

Evidence tags: **[S]** seen in this repository or in a search-result summary (the papers themselves could not be opened: arXiv is unreachable from this environment, so every literature number is a snippet, never a table read from the paper); **[C]** computed in this session (scripts named in "Evidence notes"); **[R]** recalled, unconfirmed; **[I]** inferred, with the step shown. A figure keeps the tag of its weakest input. Times shown to Boss are US Eastern (EDT) only.

## Summary

- **The rule changes in one place.** Boss's 5:39 PM EDT rule said scatters carry no systems, planets or moons. His 5:51 PM EDT message narrows that: a scatter row still carries **no planets, moons or belts**, but it **does carry its star-like companions** (stars, brown dwarfs, white dwarfs, neutron stars, black holes, pre-main-sequence stars), generated when the scatter is drawn. Planets, moons and belts stay fill-only.
- **A scatter row becomes a system key.** Its point is the system's barycentre, its mass is the system's total mass, and its members and their orbits sit in two new child tables with one column per component and no JSON (section 5). A fill builds the system from those rows instead of drawing the companions again, so a companion that landed in a scatter is the companion the sector gets.
- **The cost is small.** At the default galaxy the scatter holds about 1.13e6 star rows and 5.9e7 black-hole rows [C]. Companions add about 4.9e6 member rows (the placed object included) and 3.0e6 orbit nodes, roughly 0.6 GB against the 8 GB the base rows already need, and about half an hour of CPU at scatter time (section 6).
- **The census says what the probabilities are** (section 4). Stars follow the existing O to M table: a scattered O or B star has a mean of 3.1 or 2.3 stars (the rows' mix gives 2.77 and 2.28) [S, C]. A black hole has a companion at about 1.25% (star 1.0%, black hole 0.2%, neutron star 0.05%), not the 10% to 30% the earlier note took from three Gaia systems [R, I]. A neutron star has a neutron-star companion at about 0.04% [C, I]. Every default is one number in one tuning dict, with a range.
- **Neutron-star pairs would be missed by the mass cut, so promote the pairs.** A neutron star never passes the 14 solar mass cut, so a double neutron star would be a fill-time object. A pair is rare (about 1.7e5 of them in the default galaxy) [C, I], so the scatter places every bound compact pair (neutron star with neutron star, with black hole, black hole with black hole) whatever its mass, about 3e5 extra rows. Default on; Boss decides (section 4.5).
- **The theoretical limit is high; the observed one is 7 confirmed and 12 proposed.** A hierarchy can in principle hold a dozen levels between a contact binary and the galaxy's tidal limit [C]; seven stars is the most confirmed, and ten candidates reach 12 [S]. **The default cap is now 12** (section 16), with the Galaxy holding about 1.1e6 eight-star and 1e2 twelve-star systems [C]. Young non-hierarchical groups (trapezia) decay in 1e4 to 4e4 years [S]; they are now stored while younger than their decay age, never after it. Compact members add three filters: survival of each collapse, the gravitational-wave merger time, and the relativistic quench of Kozai-Lidov cycles (section 3.3).
- **Motion needs almost nothing new** (section 7). The barycentre point moves with the system's total mass (GEN.223); a member's position is a function of time computed on read from its stored phase, epoch and period, so the nightly update writes nothing per member.
- **Section 13 is the catalogue Boss asked for at 6:08 PM EDT**: about sixty types in six families (architecture, interacting binaries and merger products, white dwarfs, compact objects, substellar companions and disks, runaways). Each has an observed rate with a tag, a theoretical rate where no observation exists, where it occurs, whether it is drawn at scatter (S) or fill (F), and the star-model change. Most types are not separate generators: a Roche-lobe classifier turns the drawn tree into contact, Algol, post-common-envelope, accreting and merger-product pairs, and a tuning dict aligns the result to the observed fractions (13.7).
- **Twenty proposed TODO items** (section 11: eleven from the first note, nine from the catalogue) and ten rewordings of existing items. Eight decisions for Boss (section 12).

## 1. What Boss asked and how it fits the rules

| Boss's words (5:51 PM EDT) | Meaning for the design |
|---|---|
| "no forced representation in a sector (except if one landed there in a scatter before generation)" | A companion is drawn at its observed rate. A sector never gets a companion to fill a quota. A companion that a scatter placed is part of the system the fill builds. |
| "stars and star-like objects (former stars, protostars)" | A member's kind is one of: main-sequence star, giant or supergiant, white dwarf, neutron star, black hole, brown dwarf, pre-main-sequence star, protostar. White dwarfs, neutron stars and black holes are the "former stars". |
| "if a star, black hole, neutron star, etc has a companion I want it generated at scatter time" | The scatter's placed object (a bright or massive star, a black hole, a neutron star when scattered, a planetary-nebula central star, an intermediate-mass black hole) rolls for companions in the scatter pass. |
| neutron stars with neutron stars, black holes with neutron stars, other black holes, stars with neutron stars or black holes | All five pairings are in the tables of section 4, each with a source and a range. |
| "spawn a research thread to find out how many multi-star (or star-like object) stable systems are there and what the theoretical limits are and probabilities of each" | Sections 3 and 4. |
| "expand our multi-star system definitions" | Sections 4 to 8: the member kinds, the tree, the stored form, the generator. |
| "represent gravity" | Section 7: the barycentre carries the total mass; the Hill-sphere rule of GEN.109 reads it. |

The earlier decisions that still hold: scatters carry no planets, moons or belts; nebulae need stars and phenomena, never planets; planets and moons are made only when a sector is filled with everything; a fill represents every kind at its observed odds; redesign freely, the schema is not a constraint (Boss, 5:27 PM EDT).

## 2. What exists today and what a scatter row is [S, C]

- **Rows.** `bright_stars` holds a star (full star record, a `seed`, a point in mpc); `phenomenon_scatter` holds `kind`, `subtype`, a point, a `seed`, `mass_solar`, `map_luminosity_sol` (`phenomenon_scatter.py:73`, `docs/reference/database-schema.md` 1795 and 1948). A row is a point and a seed. A fill later builds the object from the seed.
- **Companions today.** A star pre-placed in `bright_stars` gets its companion at fill, when `_add_preplaced_systems` builds a `StarSystem` from the row (`run_sector.py:515`); `BinaryStarProxy` and `WideBinaryPair` draw one companion, close or wide, 0.5 each (`binary.py`, `wide_binary.py`). A scattered black hole or neutron star becomes a bare phenomenon (`run_sector.py:449`, `:454`) and `StarSystem.__init__` skips binary generation for a remnant (`system.py:255`).
- **GEN.129 and GEN.130 are unbuilt.** There is no `orbit_nodes` table or tree generator in the code (a search of `src` finds neither; they exist only in the design notes and `TODO.md`). The tree generator of `multistar-and-compact-systems.md` 5.3 is therefore the thing this note puts at scatter time.
- **How many rows the scatter holds** at the defaults (galaxy of 2.953e11 systems, stellar mass limit 14 solar masses, luminosity floor 9,000 solar luminosities, compact limit following the stellar limit) [C, from `star_population.mass_band_fraction`, `bright_band_fraction` and `phenomenon_scatter.share_above`]:

| Scatter class | Rows | Source |
|---|---|---|
| Stars born at 14 solar masses or more (mass pass) | 5.24e5 | `mass_band_fraction(14)` x systems |
| Lighter stars at 9,000 solar luminosities or more (luminosity pass) | 6.03e5 | `bright_band_fraction(9000, mass_range=(None, 14))` |
| Stellar black holes at 14 or more (40% of 5.0e-4 per star) | 5.90e7 | `share_above` 0.4 |
| Intermediate-mass black holes | 1.48e5 | rate 5.0e-7 per star |
| Planetary nebulae (white dwarf at the centre) | 2.95e4 | |
| Supernova remnants | 3.2e3 | |
| Neutron stars | 0 (all drawn at fill; they top out at 2.2 solar masses) | |

  The scatter therefore holds about 6.1e7 rows, 97% of them black holes. The star classes are small.

## 3. How many multi-star systems exist, and the limits

### 3.1 The census, as the earlier note has it and what is new

The star numbers are in `multistar-and-compact-systems.md` 3.1 to 3.3 and stand ([R3646]: 56% single, 33% double, 8% triple, 3% higher among FGK dwarfs within 25 pc [S]; [R4384] [S]; [R3108]: 0.50 companions per solar-type primary rising to 2.1 for O stars [S]; the observational ceiling is seven stars, none certain [S]). New in this note, from search-result summaries [S]:

| Quantity | Value | Source |
|---|---|---|
| Triple fraction of early B and O stars | about 50% and 70% | [R3108], as quoted in [R4203] |
| Only 16 +- 8% of the stars that dominate the core-collapse supernova rate are single | | [R3108] |
| White dwarfs: share with a companion | 25% to 26% (232 white dwarfs within 25 pc: 74% single, 26% not); the 20 pc sample gives 32% | Holberg et al.; [R4404] |
| Black-widow and redback spiders | 41 confirmed and suspected in the field; with the redbacks they are about 80% of binary field millisecond pulsars with spin period under 2 ms; companion masses 0.07 to 0.1 solar masses are almost absent | [R4279] |
| Surviving the first supernova | binary survival about 14% (+20, -10), neutron-star and black-hole remnants in similar numbers; 1% to 10% of core-collapse supernovae occur with a compact companion present | [R4441] |
| Fraction of massive-star binaries that end as a merging black-hole pair | about one in eight interacting binaries massive enough to make black holes (COMPAS, low metallicity) | [R4512] |
| Gaia dormant neutron-star candidates | 21 astrometric binaries with solar-type stars and dark companions near 1.4 solar masses; Gaia NS1 (1.90 +- 0.04 solar masses, 731 d, a halo G star) | El-Badry et al. ([R1239], 2405.00089) |
| Gaia black holes with luminous companions at AU scales | about 0.1% of Milky Way black holes (prediction of about 30 in Gaia DR4); up to about 900 detectable in the thin disc | [R2463, R3233] |
| Black-hole pairs in the Galaxy | about 1.2e6 (cosmological simulation, mean system mass 28 solar masses) | [R2539] as quoted in [R3348] |
| Double neutron stars: merger rate per Milky Way | 42 (+30, -14) per Myr ([R3557], eight merging systems), 37 (+24, -11) [R3557], 32 (+19, -9); gravitational-wave based 5.3 (+17.6, -4.9) and 11.0 (+19.2, -8.2) per Myr | [R3557], 2107.13307, 2402.04658, 2604.05059 |
| Local merger-rate density (GWTC-4.0, 90% intervals) | binary neutron star 7.6 to 250 (published 8.8 to 250), neutron star with black hole 9.1 to 84 (9.3 to 86), binary black hole 14 to 26 (13 to 26), per Gpc^3 per year | [R0006] |
| Trapezium-type systems decay in | 1e4 to 4e4 years (bound ones); after two crossing times 41 of 100 still trapezia; the ejected star is usually the lightest | [R0075], 1701.03440, 1908.11437 |
| Black-hole binaries in triples | about half of those formed in the TSE model are in triples | [R4203] |

### 3.2 The extreme members of each kind of company [R unless noted]

| System | What it shows |
|---|---|
| Nu Scorpii, probably 7 stars; AR Cassiopeiae and 65 UMa candidates | the stellar ceiling [S] |
| PSR B1913+16 (Hulse and Taylor, 7.75 h) and the double pulsar PSR J0737-3039 (2.45 h, about 0.006 AU, merger in about 85 Myr) | neutron star with neutron star exists; about twenty are known in the Galaxy [R], eight of them merge in a Hubble time [S: Pol et al.] |
| GW170817 and GW190425 | binary neutron stars merging in other galaxies [R] |
| Cygnus X-1 (a black hole of about 21 solar masses and an O9.7 supergiant of about 40, 5.6 d) | black hole with a massive star [R] |
| Gaia BH1, BH2, BH3 (periods about 186 d, 1,300 d, 4,250 d) | dormant black hole with a star [S, `multistar-and-compact-systems.md` 7.2] |
| Gaia NS1 | dormant neutron star with a star [S] |
| PSR J0337+1715 | millisecond pulsar with a white dwarf (1.6 d) and a second white dwarf (327 d): a compact triple [S] |
| V404 Cygni | black hole with a star (6.5 d) and a tertiary at about 3,700 AU: a compact triple [S, `exotic-environments...` 2.3] |
| PSR B1620-26 in M4 | pulsar and white dwarf with a planet: a captured circumbinary giant [S] |
| Neutron star with black hole | none known in the Galaxy; the gravitational-wave catalogue has several [R] |
| Black hole with black hole | none known in the Galaxy; a hundred and fifty gravitational-wave mergers [S: GWTC-4.0 has 158 events] |

### 3.3 The theoretical limits

1. **Hierarchy depth.** A level needs an outer-to-inner axis ratio of about 3 or more for equal masses (Mardling and Aarseth 2.9 to 3.3 coplanar [S, `multistar...` 5.1]) and observed period ratios are mostly above 10, i.e. axis ratios above about 4.6. The tightest stellar orbit is a contact binary at 0.01 AU and the widest bound orbit the generator keeps is 1e4 AU (20,000 AU for the root), a span of 1e6 to 2e6 in axis. A chain of axis ratio 3 fits `ln(1e6)/ln(3)` = 12.6 levels, and a ratio of 5 fits 8.6 [C]. So the stability criteria allow a chain of **a dozen or more stars**, and the earlier note's "theory allows several hundred" [S, [R1529]] needs branching. Nobody has seen more than seven. ~~Default: cap at seven~~ (`MAX_STARS`), superseded: the default cap is 12 (Boss, 2026-10-15, 6:08 PM EDT; the group cap stays 12, unified-star-like-objects.md).
2. **Non-hierarchical systems do not last.** Trapezia decay in 1e4 to 4e4 years; the lightest star is ejected and the rest harden [S]. A system older than a few times 1e5 years is therefore hierarchical or not bound. **Rule: a stored tree is hierarchical.** A trapezium-like group belongs to the protostar cores of GEN.241 (age under 1 Myr), where it may be drawn as a flagged "decaying" group later; not in this note's build.
3. **Compact members add three filters** that a star-only tree does not have:
   - *Survival of each collapse.* Mass loss unbinds a circular orbit when more than half the pre-collapse mass is lost (Blaauw) and a neutron star's natal kick (about 265 km/s) unbinds almost every wide orbit. A Monte Carlo of a black-hole progenitor (25 to 100 solar masses from the IMF) with a Sana period law and a thermal eccentricity, **only the 17% of binaries that stay out of mass transfer**, gives per massive primary [C, `mc2.py`]:

     | Black holes collapse directly above | Companion stays bound as a star | black hole pair | black hole with neutron star |
     |---|---|---|---|
     | never (momentum-conserving kick) | 0.15% | 0.26% | 0 |
     | 60 solar masses | 0.28% | 0.98% | 0 |
     | 40 solar masses | 0.70% | 2.4% | 0.001% |
     | 25 solar masses | 2.2% | 5.4% | 0.002% |

     The result moves by 20 times with the one parameter nobody knows (where direct collapse starts), but its middle matches the censuses of 3.1 to within the same factor (section 4.4). It also shows that a neutron star's kick leaves **no wide black-hole-with-neutron-star pair**: those are made only in the interacting channel (tight orbits), which is why they sit at hundredths of an AU in the observed gravitational-wave population.
   - *Gravitational-wave merger time (Peters).* A bound compact pair must not already have merged. Circular orbits that merge in 13.8 Gyr: 0.022 AU (2 neutron stars), 0.046 AU (neutron star and 8 solar mass hole), 0.096 AU (10 + 10); with e = 0.9, 0.094, 0.197 and 0.41 AU [C]. A log-flat draw between 0.003 and 100 AU puts 19%, 26% and 33% of pairs inside these axes. **Rule: reject a drawn compact pair whose Peters time at its axis and eccentricity is shorter than the system's age** (a one-line test; already decided in GEN.130).
   - *Relativistic quench of Kozai-Lidov.* As GEN.129's extension already says: a Kozai screen rejects (or marks `kl_active`) an inner pair only when the mutual inclination is 39.2 to 140.8 degrees, the Kozai time is under the age and the Kozai time is under the inner orbit's relativistic precession period.
4. **Compact members per system.** Each collapse in a wide tree is survived with probability of order 0.15 [C, `mc.py`/`mc2.py`; Renzo et al. 14% [S]], so k compact members cost about 0.15^k. Two (a pair) are the observed maximum in a pair; the compact triples known (PSR J0337+1715, V404 Cygni) hold one or two compact members. **Default: at most three compact members per tree**, drawn with the product rule; a tree of four (two compact pairs) has a probability under 1e-6 of a black-hole row and is not special-cased.

## 4. Probabilities

All defaults live in one dict, `MEMBER_STATE_BY_PLACED` (section 6), each with its range. The audit of section 9 counts every cell.

### 4.1 Member kinds and how a star-like companion gets its state

A member is drawn with an initial mass, then its state follows from the system's age with the same evolution code a single star uses (`physics/stellar_evolution.evolve_star`). A member lighter than the placed object (the usual case: the placed object is the heaviest) lives longer than it does, so it is alive: main sequence, giant, or **pre-main-sequence** if it is under 8 solar masses and its age is under its contraction time (`physics/pre_main_sequence.py`, GEN.242: derived from the energy a convective then radiative star releases, with its Hayashi temperature, birthline radius and energy factor `PRE_MAIN_SEQUENCE_ENERGY_SOURCE_FACTOR` fitted to the BHAC15 tracks of [R0254] (GEN.529); 37 Myr for a solar mass, 340 Myr for 0.2 and 0.73 Gyr for 0.1, with no 10 Myr cap, so a red dwarf still shrinking at 100 Myr keeps its pre-main-sequence class; it replaces the older `31 Myr x m^-2.5` law, which gave 10 Gyr at 0.1 solar masses; this is what makes the companions of a 5 Myr-old O star young), or a **protostar** if the system is under 0.5 Myr (only inside the cores of GEN.241; a scatter row is never that young). A **brown dwarf** companion has the desert rule (under 1% inside 5 AU, a few per cent wide; GEN.252). A **former star** companion appears only through the explicit compact-companion probabilities below, because a companion that is heavier than the placed object at birth contradicts "the placed object is the heaviest"; the compact probabilities stand for the mass-transfer reversal and the heavier-companion cases together. A remnant is allowed only if its progenitor's lifetime is shorter than the system age and the placed object's lifetime is longer (for a living 14 to 25 solar mass star the window is its last third); the table value is the mean over all ages and the generator divides it by the window's share and applies it inside the window.

### 4.2 Placed star (the mass pass and the luminosity pass)

The number of stars follows `multistar-and-compact-systems.md` 3.3, chosen by the placed star's **initial mass** class (O 16 and up, B 2.1 to 16, A 1.4 to 2.1, and so on). The rows' mix [C]:

| Pass | Rows | Class mix | Mean stars per system | Companion rows added |
|---|---|---|---|---|
| Mass pass (14 solar masses and up) | 5.24e5 | O 61%, B 39% | 2.77 | 9.3e5 |
| Luminosity pass (9,000 solar luminosities, lighter than 14) | 6.03e5 | B 100% (supergiants, bright giants and subgiants of 2.1 to 14 solar masses) | 2.28 | 7.7e5 |

Companion states (per system, on top of the tree's star count):

| Placed object | Compact companion | Value and range | Basis |
|---|---|---|---|
| Living star of 14 solar masses or more | neutron star 0.2%, black hole 0.1% | 0.3% (0.05% to 2%) | about 130 known high-mass X-ray binaries [R] against 5e5 such stars in the model, times a factor of 5 to 10 for quiet ones [I]; wide-channel survival of 0.15 to 2.2% of massive primaries before the first collapse [C] gives the ceiling |
| Wolf-Rayet star (GEN.243) | as above, companion is O or B in most cases | companions 50% (WN) and 90% (WC) [S, Dsilva et al.] | already in the Wolf-Rayet design |
| Bright giant of 2 to 14 solar masses | white dwarf 1%, neutron star or black hole 0.05% | (0.2% to 5%) | barium stars and symbiotic systems are all white-dwarf companions to giants [R]; wide neutron-star survival is near zero [C] |

A companion that is a brown dwarf or a low-mass star on an outer orbit comes from the tree. Planets, moons and belts are not drawn.

### 4.3 Placed black hole (5.9e7 stellar rows, 1.5e5 intermediate-mass rows)

| Companion | Probability for a stellar black hole of 14 solar masses or more | Range | Evidence |
|---|---|---|---|
| none | 98.75% | | |
| star (dormant Gaia-type, X-ray binary, wide) | 1.0% | 0.1% to 2.2% | 0.1% of black holes with a luminous companion at AU scales (Gaia prediction) [S]; wide-channel Monte Carlo 0.15% to 2.2% per massive primary [C]; the three Gaia systems are the only direct evidence [S] |
| black hole | 0.20% | 0.017% to 2.4% | tight pairs from the gravitational-wave rate: 1.1 to 2.2 per Myr per Milky Way, a stock of 8e3 to 1.7e4 pairs including wide ones [C], so 1.7e-4 to 3.3e-4 of 1e8 holes; the Lamberts simulation gives 1.2e6 pairs, 2.4% of holes [S]; Monte Carlo 0.26% to 5.4% [C]. The default is the geometric middle of the two |
| neutron star | 0.05% | 0.0075% to 0.07% | gravitational-wave rate of neutron-star black-hole mergers 0.8 to 7.4 per Myr, stock of 7.5e3 to 6.9e4 pairs [C]; the wide channel gives none [C] |

That is 1.25% with a companion (the earlier note's 10% to 30% [R] was an inference from three Gaia systems and does not match the censuses above or the black-hole count: 10% of 5.9e7 rows would be 6e6 dormant pairs for Gaia to have found only three). When a black hole has a companion, **30% also have an outer third member** (a star or a compact object) [I: about half of black-hole binaries are in triples in a population model, and observed massive stars are 50% to 70% triples, thinned by survival]; V404 Cygni is the model case. A black hole's tree obeys the Peters rule for any compact inner pair and the direct-collapse rule of `exotic-environments...` 2.3 for any companion beyond 1,000 AU.

Intermediate-mass black hole: star companion 1% (a captured star in a cluster or a satellite), nothing else (default; no source). Hypervelocity stars: no companion (the Hills mechanism ejects one star of a disrupted binary) [R]. Supernova remnants: the core-collapse core is built at fill with the young-neutron-star company of 4.4; the 3.2e3 rows carry no companion now.

### 4.4 Placed neutron star, white dwarf and the rest

Neutron stars are not scattered at the default cut, so these apply at fill and for the promoted pairs of 4.5. A neutron star's companion [multistar 7.3, unchanged except for the new pairs]:

| Class | No companion | Companion | New: compact companion |
|---|---|---|---|
| Young neutron star | 97% to 99% | 1% to 3% (Be or main-sequence star; wide orbits unlikely) | neutron star 0.04% (range 0.014% to 0.095%) |
| Recycled millisecond pulsar (10% of pulsars) | 20% to 25% | 75% to 80% (white dwarf, spider, rare neutron star) | the neutron-star share of those is part of the 0.04% |
| Any neutron star | | | black hole 0.0025% (range 7e-4% to 0.007%): the same pair as a black hole's 0.05% [C] |

The neutron-star number: 6.8e4 to 4.7e5 double neutron stars in the Galaxy (rate x 2.5 Gyr / 0.19 of orbits that merge in a Hubble time) is 1.4e-4 to 9.5e-4 of 1e9 neutron stars [C]; the default is their geometric middle (3.6e-4, rounded to 0.04%). Cross-check: 20 known double neutron stars among about 400 binary millisecond pulsars is 5%; with millisecond pulsars 0.2% of neutron stars and 78% in binaries, 5% of those gives 8e-5, and a second neutron star in the pair that is not itself recycled raises it to the same order [I].

White dwarf: companion 30% (22% wide), per `multistar-and-compact-systems.md` 7.3; the new cross-check is Holberg's 26% [S]. Planetary-nebula central star: close companion 20% [R]. Neither is scattered except the planetary-nebula star (2.95e4 rows).

### 4.5 Promote bound compact pairs to the scatter

Boss named neutron-star pairs. At the default cut they would never be scattered, because 2.8 solar masses is under the 14 solar mass cut: a pair's tidal sphere `1.44 pc x (2.8)^(1/3)` is 2.0 pc, inside one 4 pc sector (`phenomenon-scatter-mass-cut.md`: the cut falls where the sphere equals a sector edge), so the cut's logic says a fill covers them. That logic does not hold for what Boss wants (to see and navigate the compact binaries, to mark mergers, to have pulsars for timing), and the rows are cheap because the pairs are rare:

| Pair | Probability per member | Pairs in the default galaxy |
|---|---|---|
| neutron star with neutron star | 3.6e-4 (1.4e-4 to 9.5e-4) of 9.5e8 | 1.7e5 (0.7e5 to 4.5e5) |
| neutron star with black hole (below the cut) | 3e-4 of 8.9e7 sub-cut holes | 2.7e4 |
| black hole with black hole (both below the cut) | 0.2% of 8.9e7 | 9e4 |

**Default: the scatter places every bound compact pair whatever its mass** (about 3e5 rows, 0.5% of the black-hole rows). The below-cut fill draws then leave those objects out (section 6.4). Pairs with only one compact member (a neutron star with a star, a pulsar with a white dwarf) are not promoted: their mass is the neutron star's.

### 4.6 Orbits per pairing

| Pairing | Axis law | Eccentricity | Checks |
|---|---|---|---|
| Star with star | the truncated log-normal in period of `multistar` 4 (the O and B mean toward log P(d) = 3.5, Sana's `(log P)^-0.55` below that) | circular under 7 d, thermal above, capped at `1 - (P/2 d)^(-2/3)` | Mardling-Aarseth, Eggleton-Kiseleva or the Equation 4 switch; Kozai screen |
| Black hole or neutron star with a star, interacted (tight) | log-uniform in period from 0.3 to 1e4 days [R: Gaia BH1 to BH3 are 186 to 4,250 d, X-ray binaries hours to days] | circular under 7 d, else thermal | Roche-lobe floor; the star's radius |
| Compact pair, tight | log-uniform in axis, 0.003 to 100 AU | circular under 0.1 AU, else thermal | **Peters time against the system age**; relativistic Kozai quench |
| Any wide companion (over 1,000 AU) of a collapsed object | as stars | as stars | **only for a direct-collapse black hole** |
| Tertiary of a compact binary | as the star law, outer axis at least the Mardling-Aarseth minimum | as stars | the tree checks; V404 Cygni (axis ratio 25,700) passes |

As built (GEN.261, pinned by GEN.266's tests in `test_gen_compact_orbit_laws.py`): every row of this table is in `multistar._Builder` (`law`, `period_days`, `eccentricity`, `place`): a pairing of a compact object with a star uses the `interacted` law (`INTERACTED_PERIOD_DAYS`), a pairing of two compact objects the `compact-pair` law (`COMPACT_PAIR_AXIS_AU`, `COMPACT_PAIR_CIRCULAR_BELOW_AU`, the Peters filter), a wide companion of a hole that did not collapse directly is capped at `WIDE_COMPACT_AXIS_AU`, and every other pair, a tertiary included, takes the star law. The relativistic quench is in `stability.kozai_verdict`.

## 5. The stored form

The placed object stays the scatter row (nothing about a row without companions changes). A row with companions gains columns and child rows:

**On `bright_stars` and `phenomenon_scatter`:**

| Column | Type | Null | Notes |
|---|---|---|---|
| `member_count` | SMALLINT UNSIGNED | nullable | Number of members including the placed object; NULL for a single object. |
| `system_mass_solar` | DOUBLE | nullable | Sum of the member masses; NULL for a single object. The barycentre's mass. A reader uses `COALESCE(system_mass_solar, mass_solar)`. |

The row's point, velocity, epoch and next-update columns (GEN.223's point structure) are the **barycentre**.

**`scatter_members`** (one row per member, the placed object included as member 0 so the tree is uniform):

| Column | Type | Notes |
|---|---|---|
| `source` | TINYINT | 1 = `bright_stars`, 2 = `phenomenon_scatter`. |
| `source_id` | BIGINT UNSIGNED | The row's id. |
| `member_no` | SMALLINT | 0 is the placed object. Primary key `(source, source_id, member_no)`. |
| `parent_node_no` | SMALLINT | The orbit node this member is a child of. |
| `slot` | SMALLINT | 0 and 1 for a pair (`a`, `b`); 0 to N - 1 for a member of a `group` node (section 14). |
| `kind` | VARCHAR(20), CHECK | `main-sequence`, `giant`, `white-dwarf`, `neutron-star`, `black-hole`, `brown-dwarf`, `pre-main-sequence`, `protostar`, plus the catalogue kinds `stripped-helium`, `thorne-zytkow`, `wolf-rayet` (section 14). |
| `state`, `flags` | VARCHAR(16), SMALLINT | Added by the catalogue: the member's evolutionary state and its flag bits (section 14). Nullable. As built, `state` is the Yerkes class of a star (`V`, `III`, `VII`, ...), `PMS` for a pre-main-sequence star (the column is 16 characters wide), `protostar`, and NULL for a remnant or a brown dwarf. |
| `offset_*_au`, `velocity_*_kms`, `decay_age_yr` | DOUBLE | Nullable; only the members of a `group` node (a young non-hierarchical system). |
| `mass_solar` | DOUBLE | Current mass. |
| `initial_mass_solar` | DOUBLE | Birth mass (the progenitor's for a remnant); a star's current state is rebuilt from it and the age. |
| `seed` | BIGINT UNSIGNED | The member's own build seed (details, planets later). |

**`scatter_nodes`** (N - 1 rows for N members), the same names and units as GEN.129's `orbit_nodes` so a fill copies column for column:

| Column | Type | Notes |
|---|---|---|
| `source`, `source_id`, `node_no` | | Primary key; `node_no` 1 is the root. |
| `parent_node_no`, `slot` | | NULL parent for the root. Children are found by their parent and slot (`slot` is SMALLINT, as on members). |
| `mass_solar` | DOUBLE | The subtree's mass. |
| `separation_km`, `eccentricity`, `inclination_deg`, `ascending_node_deg`, `arg_periapsis_deg` | DOUBLE | The pair's orbit (the relative orbit of child b about child a). |
| `mean_anomaly_deg`, `epoch_unix`, `period_years` | DOUBLE | The phase at an epoch. |
| `kind` | VARCHAR(16) | `close`, `wide`, `compact` or `group` (for naming and the checks). |
| `interaction` | VARCHAR(20) | Added by the catalogue: `detached`, `semi-detached`, `contact`, `accreting`, `post-common-envelope`, `merger-product` (section 14). |
| `flags` | TINYINT | bit 0 interacted (drawn from the tight-channel law), bit 1 `kl_active`. |

One column per component, no JSON. There is no foreign key (a row table of 6e7 entries with a polymorphic parent), so the build deletes the child rows of a built row in the same transaction, and a nightly check (extends `check-db`) counts child rows without a parent.

`phenomenon_scatter_classes` gains no column: a compact pair is counted under its placed object's class, with a `members` count in the Phenomena table.

## 6. The generator and its cost

### 6.1 One function, two call sites

`draw_member_tree(placed, age, rng)` returns the members and nodes. It is GEN.129's hierarchy generator (design 5.1 to 5.3 of the earlier note: draw N from the O to M table, masses, shape, root separation, then each pair top-down under the stability criteria and the Kozai screen) with three additions: the **member kinds** and the compact-companion probabilities of section 4, the **compact filters** of 3.3 (Peters, quenched Kozai, wide-companion rule), and the **state of each member** at the system's age (4.1). The scatter pass calls it for every placed star and black hole, and the fill calls it for every object the fill creates (stage D's other systems, stage G's below-cut remnants), so the two places draw from one set of laws.

### 6.2 Streams and determinism

The tree is drawn from `draw.bound(row.seed)` with the suffix `members`, in a fixed order (GEN.55, GEN.56): the existing rule that a new stage never moves an older stage's draws. The planets later come from the suffix `planets`. A re-plan with the same galaxy seed regenerates the same members; after a plan the stored members are authoritative, so a later change to the laws changes the next scatter, not the existing rows (no drift between a stored summary and a regeneration).

**As built (GEN.262).** `generation/scatter_companions.py` draws the tree from `draw.Stream(f"{seed}:members")`: first the system's age (a compact object's from `COMPACT_REMNANT_AGE_RANGE_GY`, a planetary nebula's central white dwarf's progenitor and main-sequence life, a star's from its row), then `draw_member_tree`. The address is not an input, so a row gives the same members wherever it sits. The nodes' epoch is the orbit epoch when an orbit update has run and `tuning.SCATTER_PHASE_ORIGIN_UNIX` (2000-01-01) before that, never the database's clock, so two runs of one seed store the same nodes. The row draws of the passes are untouched: the passes yield exactly the rows they yielded before, and a tree is drawn from the row at the write.

### 6.3 Bulk screening

A black-hole scatter row is made in bulk by numpy (`_band_rows`, about 1 microsecond a row). The companion decision must not turn that into a Python loop over 6e7 rows. Therefore the bulk pass draws **one uniform per row** from the shell's stream and compares it with the cumulative companion probability of the row's class (1.25%); only the hits (about 7.4e5 black-hole rows) call `draw_member_tree`.

**As built (GEN.262).** The uniform is not drawn from the shell's stream, because a new draw there would move every later draw of the shell (6.2's rule). It is a pure function of the row's seed, `scatter_companions.screen_unit(seed)`: the SplitMix64 finalizer over the seed and a per-question salt, integer arithmetic only, about a microsecond. A row it does not name costs that and a comparison with the class's total odds; a hit gets its companion state from the same uniform (`multistar.companion_state_at`) and calls `draw_member_tree(..., companion=state)`, so the tree takes the state the screen named. A neutron star is a millisecond pulsar when a second uniform of the seed is under `NEUTRON_STAR_PULSAR_CHANCE * PULSAR_MILLISECOND_CHANCE` (0.2%), and then takes the recycled-pulsar table. Star rows are not screened: their tree always runs, since most O and B stars have company (3.3). A planetary nebula's central star is a white dwarf of the `planetary-nebula-star` table; supernova remnants, hypervelocity stars, quasars and the nucleus carry no tree.

### 6.4 Density bookkeeping: a row is a system

The rates in `PHENOMENON_DENSITY_PC3` count **objects**. A black-hole pair would be two objects and, if both were placed independently, a double count, and a pair of which one member is under the cut would be counted once by the scatter and once by the fill. Rules:

- A scatter row is a **system keyed by its placed object**: its rate is the object rate divided by `1 + p_placed`, where `p_placed` is the expected number of **same-class** companions per row (0.20% for a black hole; 0.04% for a neutron star), so the black-hole density is unchanged. The correction is 0.2%, below the rate uncertainty, but it is a rule so the audit can check it.
  **As built (GEN.262):** `phenomenon_scatter.row_density_divisor(kind, subtype)` is `1 + p_placed` (1.002 for a stellar black hole; 1 + 0.998 x 1.8e-4 + 0.002 x 2.0e-3 for a neutron star, from the table rather than the rounded 0.04%; 1 for every class whose companions are not of its own class, a planetary nebula's white dwarf included), and `_scattered_rates` divides each class's rate by it. The below-cut fill draw is unchanged for now: it builds no companions until stage G draws trees (GEN.251).
- A star row is keyed by its heaviest member, which is what the census (a per-primary census) counts already.
- The fill's below-cut draw for a class leaves out objects that are members of a placed system: the rate of the fill draw is the object rate times `(1 - p_member_of_placed)`, about 0.04% for neutron stars. Stated for the audit; it changes nothing a user sees.
- A promoted compact pair (4.5) is a row with the placed object the heavier member.

### 6.5 Cost [C]

| Item | Count | Cost |
|---|---|---|
| Star rows with companions (O 90%, B 75%; mass pass 61% O and 39% B gives 0.84 of its rows, luminosity pass 0.75) | about 8.9e5 systems; 1.7e6 companions | N = 2 to 4 in under 1 ms, N = 7 in 53 ms (average) [S, `multistar` 5.3]; about 1.2 ms mean, so **about 20 CPU minutes** |
| Black-hole rows with companions | 7.4e5 systems; 9.6e5 companions (7.4e5 plus 30% tertiaries) | about 0.3 ms each, **about 4 CPU minutes** |
| Promoted compact pairs | 3e5 systems | **about 2 CPU minutes** |
| Member rows (companions 2.96e6 plus the 1.93e6 placed objects) | about 4.9e6 | 55 bytes each: 0.27 GB |
| Node rows (one per companion) | about 3.0e6 | 110 bytes each: 0.33 GB |
| Row columns | 6.1e7 rows | a nullable DOUBLE and SMALLINT that are NULL cost under a byte each; under 0.1 GB |
| Base rows today | 5.9e7 x 138 B + 1.1e6 star rows | about 8 GB |

So the scatter grows by about 7% in space (0.6 GB for the child tables) and by about half an hour of CPU in time. A cut that placed every neutron star (`--compact-min-mass 1`, 1.2e9 rows, 161 GB) would not be needed for neutron-star pairs because of 4.5.

## 7. Motion (GEN.223, GEN.109) [S, C]

- **The system is one object in space.** The barycentre carries the position, velocity and epoch columns of GEN.223's point; `system_mass_solar` is its mass. The nightly update moves it with the rest of the galaxy. Nothing new in the update algorithm: the Hill-sphere influence rule of GEN.109 reads `system_mass_solar` and gets a more massive body where there is a pair.
- **No per-member writes.** A member's position is `barycentre + offsets from the tree`, with each node's mean anomaly at time t computed as `mean_anomaly_deg + 360 (t - epoch_unix) / period_years` (mod 360) on read, as `positions_at` already does for a built pair. The nightly update therefore touches no `scatter_nodes` row; 3.0e6 node rows cost nothing per night. When a fill builds the system, the node rows copy to `orbit_nodes` with their epoch and phase, and the built system's own phase advance (`advance_orbital_phases`) takes over.
- **As built (GEN.265).** `planetgen/physics/system_tree.py` computes every member's offset from the barycentre at any time from the stored nodes (Kepler's equation for each node, the two children either side of the node's centre in proportion to the other's mass, a `group` node's members flying straight from the node's epoch), and `planetgen/db/scatter_systems.py` stores, reads and places a tree under a scatter row (`member_positions_pc` is the barycentre's stored point plus the offsets). `store.advance_scatter_points` is unchanged: it moves the barycentre like any scatter point and never touches `scatter_members` or `scatter_nodes`. `scatter_systems.point_masses_solar` is the mass a point carries (`system_mass_solar`, else the object's own) for the influence rule of GEN.109.
- **The position resolution hides the tree.** A stored point is in whole milliparsecs, 206 AU. The barycentre-to-primary offset of the widest tree (10,000 AU at a mass ratio 0.5) is 0.024 pc [C], 24 mpc, well inside a 4 pc sector and smaller than the mapping detail anyone looks at; a close tree is a fraction of one stored unit. The row's point is the barycentre and no member needs its own point until the system is built.
- **Natal velocities.** A compact pair's centre-of-mass velocity is the placed object's (a black hole's is the galactic velocity of its class); the kick of the collapse that made a bound pair is already in the survival rule and does not add a second velocity.
- **Kozai-Lidov and secular effects** inside the tree are ignored over game time, as `multistar` 6.4 decided; the quench rule only decides which trees are kept.

## 8. How a fill reuses the members

The stages of `phenomena-catalogue-and-fill.md` 3 stay; this changes what stage A, B and C do:

- **Stage A and B (pre-placed stars and scattered rows become objects).** The system is built from the stored members and nodes: the placed object is member 0, the others are built from their `seed` and `initial_mass_solar` and the system age, and the nodes are copied into `orbit_nodes`. `StarSystem` takes a `member_tree` argument that **skips its own companion draw** (today it draws one and skips it for a remnant, `system.py:255`). A failed build (no room) drops the outermost body and then, if need be, the outermost member, never the placed object, as the catalogue note says.
- **Stage C (the system pass) shrinks.** It no longer draws companions for a scattered object: it draws **planets, moons and belts only**, per the kind's planet rules (no planets under 10 Myr, millisecond-pulsar planets at 0.7%, Holman-Wiegert zones taken from the stored nodes), from the stream `seed + planets`. For an object that has no stored members (a below-cut remnant, a fill-created star) the same stage draws the member tree first, with the same function, then the planets.
- **Stage D (the cell's other systems) and stage G (below-cut remnants)** call `draw_member_tree` for every object they create, so a system made at fill has companions at the same odds as one made at scatter time.
- **Nothing is forced.** No stage adds a companion to meet a quota. A sector holds, among the objects the fill creates, companions at the table rates; among the pre-placed ones, whatever the scatter drew.
- **Cleanup.** A built row's `scatter_members` and `scatter_nodes` rows are deleted in the build's transaction; the built system's own tables are then the only record.

## 9. Audit and tests

Extends the representation audit of the catalogue note (TEST.160, TEST.161):

1. **Companion states.** Over a seeded scatter of at least 1e6 rows, the share of rows with each companion state (none, star, black hole, neutron star, white dwarf) lies within four Poisson sigma of the table; the mean stars per system of each class match `multistar` 3.3 within two sigma.
2. **Row bookkeeping.** Object density per class before and after members matches (the `1 + p_placed` rule); no object is both a row and a member.
3. **Stability.** Every stored tree passes Mardling-Aarseth, Eggleton-Kiseleva or Equation 4 (the switch), the Kozai screen and the Peters rule.
4. **Determinism.** Two scatters with the same seed store identical members; the planets of a built system do not change when the member laws change for the next scatter.
5. **No planets in a scatter.** Property test: nothing in `scatter_members` or `scatter_nodes` is a planet, moon or belt.
6. **Fill reuse.** The stars of a built system equal the stored members; the nodes equal the stored nodes; the child rows are gone after the build.
7. **Motion.** After a nightly update the barycentre point has moved and no `scatter_nodes` row has been written; a member's position at time t equals its position computed from the stored phase; a round trip of the point columns is exact (GEN.223's own test).
8. **Examples as unit tests.** The double pulsar (2.45 h, 1.34 + 1.25 solar masses) passes Peters (85 Myr is below a 10 Gyr age, so a drawn pair with that axis is *rejected* for an old system and accepted for a young one); V404 Cygni (9.7, 6.5 d and 10.6 at 3,731 AU) passes the stability test; a wide companion of 3,000 AU to a 14 solar mass hole that is not direct-collapse is rejected.
9. **Census anchor test.** Over the full default scatter the counts of stellar-mass black holes with a black-hole companion lie between the gravitational-wave stock and the Lamberts count (1.7e-4 to 2.4e-2 of holes).

## 10. What this changes in the earlier notes and TODO

| Where | Was | Now |
|---|---|---|
| `GEN.129` last paragraph ("multi-star systems around a scattered object appear at sector fill only") and `GEN.130` last paragraph | Remnant systems only at fill, by promotion through the scatter | Companions at scatter; planets and moons at fill |
| `GEN.251` system pass | Rolls one system (companions, planets, belts, moons, disk) for every scattered object at fill | Rolls planets, moons, belts and disks only for scattered objects; companions come from the stored members. Objects made at fill still draw both |
| `phenomena-catalogue-and-fill.md` section 5 | Stellar black hole `P(system)` 0.20 | Companion 1.25% (star 1.0%, black hole 0.2%, neutron star 0.05%) |
| `phenomena-catalogue-and-fill.md` section 3 rule "Never planets for a scatter row alone" | A scatter row stores a point and a seed | A scatter row stores a point, a seed and its members; nothing with a planet, moon or belt |
| `multistar-and-compact-systems.md` 7.3, black-hole row | "10 to 30% [R]" bound stellar companion | 1.25% in total (this note, 4.3); the other rows unchanged |
| `GEN.223` | Scatter placeholders get point columns | The point is the barycentre of a system, with `system_mass_solar`; members need no nightly write |
| `GEN.109` | Mass of the body in the point-mass table | `COALESCE(system_mass_solar, mass_solar)` |
| `exotic-environments...` Peters rule and wide-companion rule | In GEN.130 slice (a) | Moves into the generator of 6.1 (used by both call sites) |

## 11. Proposed TODO items

The TODO thread numbers them and assigns the lanes. Items 1 to 11 are already filed (merged #1651): item 1 is GEN.260, 2 is GEN.261, 3 is GEN.262, 4 is GEN.263, 5 is GEN.264, 6 is GEN.265, 7 is GEN.266, 11 is GEN.267, 8 is UX.114, 9 is TEST.162 and 10 is DOC.50. Items 12 to 20 and the changed existing items are the catalogue's additions and wait for the TODO thread. Order is build order.

1. **Scatter member columns and tables** (Core; with or immediately after the shared GEN.223 and GEN.230 migration): `member_count`, `system_mass_solar` on `bright_stars` and `phenomenon_scatter`, the `scatter_members` and `scatter_nodes` tables of section 5, the cleanup check in `check-db`. Prerequisites: GEN.223.
2. **The member-tree generator** (Core): `draw_member_tree` of section 6.1, i.e. GEN.129 steps (a) to (d) plus the member kinds, the companion-state dict `MEMBER_STATE_BY_PLACED`, pre-main-sequence and brown-dwarf states, and the compact filters (Peters, relativistic Kozai quench, direct-collapse wide-companion rule). Absorbs the code half of GEN.130 slice (a). Prerequisites: GEN.250.
3. **Companions in the scatter passes** (Core): the mass pass, luminosity pass and phenomenon pass draw members (stars, black holes, intermediate-mass black holes, planetary-nebula central stars), with the bulk screening of 6.3 and the `1 + p_placed` density bookkeeping of 6.4. Prerequisites: items 1 and 2.
4. **Compact-pair promotion** (Core; decision 3 below): scatter every bound compact pair whatever its mass, and remove its members from the fill's below-cut draws. Prerequisites: item 3.
5. **Fill builds from stored members** (Core): `StarSystem(member_tree=...)`, stage C reduced to planets, moons and belts, stage D and G call the generator, build deletes the child rows. Rescopes GEN.251. Prerequisites: items 1 to 3, GEN.251.
6. **Nightly motion of a system** (Core, with GEN.223): barycentre point with `system_mass_solar`; analytic member phase on read; no per-member writes; GEN.109 reads the system mass. Prerequisites: GEN.223, item 1.
7. **Compact-binary rules** (Core; replaces GEN.130 slices (a) and (b) in part): the orbit laws of 4.6 for black hole or neutron star with a star, tight compact pairs, tertiaries; pulsar-planet presets stay in GEN.130 slice (b). Prerequisites: item 2.
8. **Scatter-row pages show members** (Interface): the sector and bright-star pages list the stored members of an unbuilt row and the map shows a system mass; a `Members` column in the Phenomena table. Prerequisites: item 1.
9. **Multiplicity audit** (Quality): section 9 tests 1 to 9 in the test suite, with the admin report of the representation audit. Prerequisites: items 3 and 5.
10. **Docs** (Quality): `database-schema.md` (the two new tables, the two new columns), `multistar-and-compact-systems.md` 7.3 corrected, `phenomena-catalogue-and-fill.md` sections 3 and 5 corrected, a pointer from `interstellar-object-rates.md`.
11. **Cap and directive** (Core, small): `MAX_STARS` setting (default **12**, raised by the catalogue; section 16) and a `--stars N` directive, so a user may force N or raise the cap; a prevalence knob for the N >= 5 columns, as GEN.129 already has.

    **As built (GEN.267), the settings and their plumbing; the host that reads them is GEN.630's.** `SystemConfig.STARS` (None, or the forced count) and `SystemConfig.MAX_STARS` (`tuning.MAX_STARS_DEFAULT`, 12) travel in the config's JSON, the run's plan arguments and the `system_configs` row (migration 0134, schema v134: `stars`, `max_stars`, `prevalence_multi_star`), so a stored config regenerates with the same count. `--stars N` is for one system (`planetgen system`, the system file's `stars`; sector and galaxy runs refuse it), `--max-stars N` sets the cap for any run (the system file's `max_stars`, the Generate page's Most stars in a system), and `multi_star` is a prevalence feature (`--prevalence multi_star=PERCENT`, a Generate page field that starts at the usual share): `prevalence.multi_star_scale(config)` is the factor on the `star-count` table's five-or-more columns and on the 8-to-12 tail, which the host passes as `draw_member_tree(..., max_stars=config.MAX_STARS, n_stars=config.STARS, prevalence=multi_star_scale(config))`; -100 gives none. The refusals are `run_system.stars_problem` and `max_stars_problem`, in the command line and on a system file: a count under 1, a count over the cap (the message names `--max-stars`), a cap over `tuning.MAX_STARS_CEILING` (64, a stand-in: see below), `--stars 1` with `+binary_system` or `+wide_binary`, and two or more stars with `-binary_system`. `--stars 3` with `+binary_system` stands: the system has a companion, and what the rest of the tree is belongs to the host. The usual share of five or more stars, `tuning.PREVALENCE_DRAW_SHARES["multi_star"]` = 0.46% of systems (0.44% before the mass ratios of GEN.951 moved the primary mass function), is the table's columns averaged over the primary mass function (GEN.928), because no run builds five stars until GEN.630 lands; remeasure it from generated systems then (a test holds the figure against the table until it is). The ceiling of 64 is a stand-in: a chain of axis ratio 3 fits a dozen levels in the span of axes (section 3 item 1), and branching could hold a few dozen. `StarSystem` does not read `STARS`, `MAX_STARS` or the prevalence yet, so a forced count changes nothing in a generated system until the host reads it; no generator epoch, because the defaults reproduce today's draws.

**Added by the catalogue (6:08 PM EDT)**, in build order after item 2:

12. **System-type registry and rates dict** (Core, with GEN.250): `SYSTEM_TYPE_RATES` of section 14, one entry per row of 13.1 to 13.6, plus the completeness test (every row has a generator and an audit count). Prerequisites: GEN.250.
13. **Interaction states in the star model** (Core): the Roche-lobe classifier, the `interaction` and member `state` and `flags` columns (migration folded into item 1), `mass_reversed`, `mass_gainer`, rejuvenated age, `stripped-helium`, `thorne-zytkow` and `wolf-rayet` kinds. Prerequisites: items 1 and 2.
14. **Close-binary types** (Core): contact, Algol, spectroscopic and eclipsing flags, the tertiary share of contact binaries, massive over-contact binaries, merger products (massive, field, double white dwarf bump), Be stars with stripped companions, the tight-tail weight of 13.7. Prerequisites: item 13.
15. **White-dwarf binary types** (Core; extends GEN.253): close double white dwarf, post-common-envelope pair, cataclysmic variable, AM CVn, symbiotic and barium stars, white dwarf with brown dwarf. Prerequisites: item 13, GEN.253.
16. **Accreting compact binaries** (Core; GEN.130 slice c): high-mass and low-mass X-ray binaries, Be/X-ray, microquasar flag, accretion luminosity and the Eddington limit, spider millisecond pulsars (black widow, redback) replacing the old 10% "other". Prerequisites: items 7 and 13.
17. **Young groups and trapezia** (Core, with GEN.241 and GEN.244): the `group` node, offset and velocity columns, decay age, the decayed-remnant draw, and the sector-page lifetime note. Prerequisites: items 1 and 2.
18. **Cluster-hosted compact binaries and runaway bookkeeping** (Core; GEN.159 to GEN.164): cluster rates for black holes with stars and millisecond pulsars, `runaway_origin`, the compact companion left behind by a supernova runaway. Prerequisites: item 12.
19. **Substellar and circumbinary types** (Core; GEN.252, GEN.94, GEN.95): wide and close brown-dwarf companions, brown-dwarf pairs, planetary-mass companions, circumbinary disk descriptor and planet zone. Prerequisites: item 12.
20. **Type audit and docs** (Quality): section 15 tests 10 to 14, the catalogue table in `database-schema.md` and `multistar-and-compact-systems.md`, the `theoretical` flag shown on the sector page. Prerequisites: items 12 to 19.

**Changed existing items**: item 1 (the migration also adds `state`, `flags`, `interaction`, SMALLINT `slot`, the `group` offset and velocity columns); item 2 (the generator also produces the catalogue's types); item 11 (cap 12); GEN.129 ("up to seven" becomes "up to twelve" and the trapezia sentence); GEN.130 slice (c); GEN.242, GEN.243, GEN.252, GEN.253 (each gains a pointer to its catalogue rows); `phenomena-catalogue-and-fill.md` kind list (new kinds).

## 12. Decisions for Boss, with the defaults taken

1. **Black-hole companion fraction 1.25%, not 10% to 30%** (4.3). The old number was an inference from three Gaia systems; the censuses support 1.25% (range 0.2% to 7%). Default: 1.25%. The tuning dict allows any value.
2. **The cap is 12** (3.3, 13.1, 16). ~~Seven stays the cap~~ [Edit 6:08 PM EDT catalogue]: ten candidate systems reach 12 and theory allows more. Default: cap 12, directive to raise it.
3. **Promote every bound compact pair to the scatter** (4.5): about 3e5 rows, makes neutron-star pairs exist at scatter time. If off, neutron-star pairs appear only when their sector fills. Default: on.
4. **Young non-hierarchical groups are stored; old ones never** (3.3, 13.1, 16). ~~Hierarchical trees only~~ [Edit catalogue]: a trapezium is drawn while younger than its decay age (1e4 to 4e4 yr), at fill in star-forming cores and rarely at scatter; an older group is replaced by its remnant. Default: on.
5. **Delete a row's members and nodes when its sector builds it** (8) rather than keeping two copies. Default: delete.
6. **Pair rates from the middle of the evidence** (4.3, 4.4): the black-hole pair rate is the geometric middle of the gravitational-wave stock and the Lamberts simulation, which differ by 80 times. If you want the lower or the higher end, it is one number. Default: middle (0.2%).
7. **Emergent types, tuned to observation** (13, 13.7): contact, Algol, post-common-envelope and accreting pairs come out of the Roche classifier, with one tuning knob per family so the audit counts match the observed fractions. Default: emergent plus tuning, not hand-set fractions.
8. **Theoretical-only types are included and flagged** (13.2, 13.4, 16): Thorne-Zytkow objects (100 in the Galaxy), magnetar binaries (0), quadruples with two compact pairs, planetary-mass companions. Default: generated at the predicted rate with `theoretical` shown on the page; set the rate to zero to remove one.

## 13. Catalogue of system types (Boss, 6:08 PM EDT: "every observed and theoretical type")

Boss approved the plan and widened it: plan for every observed and theoretical type of multi-star or star-like-object system, build each into the generator, and match observed probabilities as closely as possible ("observed are best, theoretical are good too"). This section is that catalogue. The earlier sections stay the mechanism (a tree of pairs, a conditional-fraction dict, one generator for scatter and fill); this section lists what the mechanism must be able to produce.

**How to read a row.** *Observed* is the real rate with its tag and the **default** the generator takes, stated against a reference class (the set of objects the fraction applies to). *Theory* is the predicted rate where no observation exists, with a confidence (high, medium, low). *Where* is the population or region. *Generated and stored* says whether the type is drawn at scatter time (S) or only at fill (F) and which stored form it uses. *Star model* is what changes in the star model. Galaxy counts use the default galaxy of 2.953e11 systems [C]. Every default is an entry in the dict `SYSTEM_TYPE_RATES` (section 14), changeable without code.

**Two mechanisms, not forty.** Most types are not separate generators:
1. **Emergent from geometry.** After the tree draws a pair's axis, eccentricity and masses, a classifier computes each member's Roche-lobe filling factor at the system's age from the stellar evolution model and assigns the pair an `interaction` state: `detached`, `semi-detached` (Algol), `contact`, `accreting` (a white dwarf, neutron star or black hole receives mass: cataclysmic variable, AM CVn, X-ray binary), `post-common-envelope`, or `merger-product`. A pair's state is a function of its birth period and the evolution, so contact binaries, Algols and cataclysmic variables appear at the rate the period law and the evolution give, and the audit checks that rate against the observed one. Where the emergent rate falls short (it will, because the generator has no common envelope or magnetic braking), the **tight-tail weight** of the period law is raised until it matches (below, "tuning").
2. **Conditional fractions.** Types that need a history the generator lacks (a post-common-envelope pair, a stripped star, a double white dwarf, a compact pair) are drawn as conditional fractions of their placed-object class (section 4), with the axis from a stated law and the Peters, Kozai and survival filters.

### 13.1 Architecture of stellar systems

| Type | Observed (default) | Theory and confidence | Where | Generated and stored | Star model |
|---|---|---|---|---|---|
| Spectroscopic / close detached binary, a under 0.26 AU | 6.6% of binaries [C, log-normal period law] (`multistar` 4) | high (period law is [R]) | everywhere; more in the thick disk and halo (close fraction 40% to 53% at [Fe/H] -1 to -3 [R]) | S and F, one node, kind `close`. Eclipsing is derived: the pair eclipses when the inclination is within `asin((R1+R2)/a)` of edge-on | none |
| Intermediate binary, 0.26 to 50 AU | 43.3% of binaries [C] | high | everywhere | S and F, node kind `close` or `wide` by the 50 AU rule | none |
| Wide binary, 50 to 10,000 AU | 43.6% of binaries [C] (GEN.448: the wide range now runs from 0.25 AU to the derived limit; 43.6% is the share of the old 50 to 10,000 AU range) | high | disk; thinned in the bulge, clusters, halo (`multistar` 6.6) | S and F, kind `wide`; planets of each star at fill | none |
| Ultra-wide pair, 10,000 to 20,000 AU (0.05 to 0.1 pc) | 6.5% of binaries [C]; Gaia finds pairs to about 1 pc [R] | medium: unbound by passing stars in 1 to 5 Gyr [R] | disk only; fraction falls with density | S and F, kind `wide`, flag `fragile`; the cap is 20,000 AU for the root, a directive raises it | none |
| Triple 2+1 (a close pair with a wide third) | 8% of FGK, 22% of B, 30% of O (`multistar` 3.3); 50% and 70% for B and O in a recent summary [S] | high | everywhere; more in the young disk | S and F, two nodes | none |
| Quadruple 3+1 (a triple with a wide fourth) | 60% of quadruples [I] | medium (shape split is a design choice) | everywhere | S and F | none |
| Quadruple 2+2 (two close pairs bound wide) | 40% of quadruples [I]; 4% of F and G systems [S, [R4384]] | high for the 4% | everywhere | S and F | none |
| Quintuple to septuple | N = 5, 6, 7: 2.4e-3, 2.9e-4, 3.0e-5 of systems [S, C]; seven-star systems confirmed: none; Nu Sco probable [S]; a septuple protostellar system in NGC 6334IN found in 2025 [S] | high that they exist | everywhere; richest in young massive regions | S and F | none |
| Eight to twelve stars | none confirmed; ten candidates with multiplicity up to 12 [S, [R1529]]; theory allows several hundred [S] | medium (stability allows 12 or more levels, 3.3) | young and massive regions | S and F, same recursion. Tail: ratio 0.12 for N = 8 and 0.1 per step after: 1.1e6, 1.1e5, 1.1e4, 1.1e3 and 1.1e2 systems for N = 8 to 12 in the Galaxy [C]. **Cap 12** (decision 2) | none |
| Young non-hierarchical group (trapezium) | 412 catalogued trapezium-type multiples, about 10% true after projection, 3 to 5 per OB association; Orion's Trapezium has 4 to 5 members; most are expanding; lifetime 1e4 to 4e4 yr bound, 41 of 100 remain after two crossing times [S]. **Default: 5% of massive (O and early B) systems younger than 0.3 Myr** [I] | high that they exist and decay | star-forming cores, OB associations (GEN.241, GEN.244) | **F** in cores; **S** only for the few mass-pass stars younger than 0.3 Myr. Stored as a `group` root node whose members carry offset and velocity columns (section 5); no Kepler elements. A group older than its decay time is never drawn; see the next row | members flagged `unstable`; the sector page shows the group with its remaining lifetime |
| Decayed trapezium remnant | the lightest member is ejected, the rest harden into a hierarchy [S] | high | clusters, associations, runaway stars | the generator draws a group older than its decay time as its **result**: an ejected lightest star (a runaway row, below) plus a tighter hierarchical remnant | the ejected star carries `runaway_origin` = dynamical |
| Multiple with a circumbinary or circumstellar disk | most stars under 3 Myr have disks [S]; Perseus protostars [S] | high | young | F; descriptor on the node (GEN.242), no bodies | `disk` descriptor |

### 13.2 Interacting binaries and merger products (stellar)

| Type | Observed (default) | Theory and confidence | Where | Generated and stored | Star model |
|---|---|---|---|---|---|
| Contact binary (W UMa, over-contact) | **0.2% of FGK dwarfs** (1/500); older estimate 1/130; up to 1/30 among blue stragglers [S]. Tertiary signature in 20 of 66 (30%) [S]; resolved tertiaries 9% of ellipsoidal binaries, 29% at the shortest periods, 65% at periods under 0.3 d [S]. Default tertiary 60% | high; massive over-contact binaries exist (VFTS 352, V382 Cyg [R]) at a low rate, default 1% of O and B binaries [I, low] | disk, older | S (massive) and F. Emergent: both Roche fills at the age, tight tail of the period law, `interaction` = `contact`. A tertiary is drawn with 60% (the Kozai channel) | contact flag, tidal locking, one shared envelope temperature |
| Semi-detached (Algol) | **0.1% to 0.2% of disk stars**, a lower limit [S]. Among eclipsing binaries: 62% detached, 8% semi-detached, 23% contact, 7% uncertain (Kepler) [S]; the LMC mix is 63/25/6 [S] | high | disk, intermediate age | S and F. Emergent: donor (subgiant, lighter) fills its lobe; the *mass-reversal* flag marks a donor that was the heavier at birth, which is how the companion of a living massive star becomes the lighter | `mass_reversed` flag; donor radius from the evolution model |
| Blue straggler and merger product, field | Triple mergers explain 20% to 25% of blue stragglers; 4% of simulated triples merge [S]; clusters 10% to 30% of cores [R]. **Default: 0.5% of old and intermediate stars** [I, low] | medium | old disk, thick disk, halo, clusters (GEN.162) | F; a single star flagged `merger-product`, age rejuvenated (set to a fraction of the system age) | `merger-product`, rejuvenated age; a blue straggler in a cluster uses GEN.162 |
| Massive-star merger product and mass gainer | about 20% of massive main-sequence stars have gained mass and about 8% are mergers [R]; 70% of O stars exchange mass with a companion [S] | medium | young, OB | S: a massive row gets state `merger-product` with 8% and the tree has one star fewer; `mass_reversed` on 20% of its pairs | `merger-product`, `mass_gainer` |
| Hot subdwarf and stripped helium star | most hot subdwarfs are in binaries or merger products [R]; stripped companions of Be stars [S] | medium | old disk (sdB), young (stripped) | F; a conditional state of a post-common-envelope pair: companion mass 0.3 to 0.5 solar masses, `kind` = `stripped-helium` | new kind `stripped-helium` (hot, small) |
| Be star | about 20% of B stars [R]; the binary fraction of Be stars is lower than B stars below 0.04 arcsec and their companions are mostly stripped stars, white dwarfs or neutron stars [S] | medium | young disk | S for mass-pass stars and F; flag `be` with a decretion disk; 30% have a stripped companion, a few a neutron star (Be/X-ray, below) | `be` flag |
| Wolf-Rayet binary | WN 52%, WC 90% [S, Dsilva] | medium | young | S (the star is a mass-pass row, GEN.243) and F; companion O, B or compact | W state |
| Luminous blue variable and hypergiant | tens in the Galaxy [R] | low | young | S and F as massive stars | GEN.246 |
| R Coronae Borealis and FK Comae stars | a few hundred and a few dozen [R] | low (merger of a white dwarf with a helium star; spun-up merger of a contact binary) | old | F, single, flagged | `merger-product` subtype |
| Thorne-Zytkow object (a neutron star inside a giant) | **none confirmed**; HV 2112 and HV 11417 are candidates and HV 2112 is probably a normal AGB star [S]; 20 to 200 predicted [S, catalogue note] | low (never observed) | old | F, flagged `theoretical`; default 100 in the Galaxy (3.4e-10 of systems): a giant with a neutron-star core, no planets, no companion | new kind `thorne-zytkow` (a red supergiant with an exotic core) |

**As built (GEN.270).** Every close-binary type of this table is drawn inside `multistar.draw_member_tree`, so the scatter and GEN.631's fill read the same tree, and nothing is stored beyond the member flags and the interaction the tree already had (the new bits `MEMBER_R_CRB` 512 and `MEMBER_FK_COMAE` 1024 are in `scatter_members.flags`, a SMALLINT; the rejuvenated age is not stored, it is the system age times a fraction drawn from the member's seed).

- *Contact and Algol.* The Roche classifier of GEN.269 decides them; GEN.270 adds the gates of the table's "Where" column (`_gate_forbids`, audit 11): a contact binary of a primary under the massive limit needs an age of 0.1 Gy or more (`tuning.YOUNG_POPULATION_MAX_AGE_GY`), an Algol needs less than 3 Gy (`OLD_POPULATION_MIN_AGE_GY`), and neither holds a brown dwarf, a pre-main-sequence star or a protostar (a contracting object fills its lobe without being an interacting binary). A pair the gate forbids is redrawn like any pair that fails a screen. `mass_reversed` goes on the donor of a semi-detached pair that was the heavier at birth; the lighter star is a `mass_gainer`. The tight tail is the knob `tight_tail_weight`, now applied in `_Builder.period_days`: it multiplies the density of periods under a day for primaries under the massive limit, and was calibrated to 6.0. Measured over 300,000 seeded systems: contact in 0.49% of FGK systems (window 0.2% to 0.77%) and Algol in 0.09% of a disk of constant star formation (0.05 to 10 Gy; the catalogue's 0.1% to 0.2% is a lower limit, and a higher weight puts more detached pairs under a day, 1.0% of systems at 6.0 and 1.35% at 8.0, so the lane stopped short of it).
- *Contact tertiaries.* The stability screen leaves the inner pair of a hierarchy only short periods, so 88% of the contact binaries the draw makes already have a third star (94% before the mass ratios of GEN.951, 84% before GEN.652). `_contact_pair_alone` keeps the others of a tree whose contact pair holds the placed star with the odds `tuning.CONTACT_TERTIARY_KEEP_SHARE` (0.48, measured over 120,000 seeded systems of 0.8 to 1.2 solar masses at 5 Gy: 34% with none kept, 88% with all kept; it was 0.26 with 48% and 94% before GEN.951) and otherwise places the contact pair alone, on the orbit it was drawn on with its orientation redrawn; a forced star count is never touched. The share comes out at 60%.
- *Eclipsing and spectroscopic.* Derived, never stored: `MemberTree.eclipsing_nodes()` (a pair whose inclination is within asin((R1 + R2) / a) of edge-on, `pair_eclipses`) and `MemberTree.spectroscopic_nodes()` (a pair closer than `tuning.SPECTROSCOPIC_MAX_AXIS_AU`, 0.26 AU). About 0.9% of systems have an eclipsing pair and 4% a spectroscopic one.
- *Merger products and mass gainers.* `draw_merger_product` (knob `merger_product_fraction`): a young star of 9 solar masses or more (`MERGER_MASSIVE_MIN_MASS_SOL`, a stand-in for the O-star edge) is a merger product with 8% and its tree has one star fewer; any other star of 0.1 Gy or more is a blue straggler with 0.5%. Both carry `MEMBER_MERGER_PRODUCT` and a rejuvenated `effective_age_yr` (0.1 to 0.6 of the system age), and the state, the radius and the lifetime are those of the star at the effective age. `draw_mass_gainer` flags 20% of the young massive stars `mass_gainer`. The six uniforms of the star-type draws are taken first, in a fixed order (`STAR_TYPE_DRAWS`), so a later draw of the tree does not move with which type was chosen.
- *Be stars.* `draw_be_star`: a B star younger than 0.1 Gy is a Be star with 20% (flag `be`); 30% of them have a stripped helium companion of 0.3 to 0.5 solar masses (`stripped-helium`, born 1.2 to 2.0 times the Be star's mass, radius from [R1875] for an sdB; period 30 to 1000 days, a stand-in) and the Be star is a `mass_gainer`. The hot subdwarfs of post-common-envelope pairs are GEN.271's.
- *R Coronae Borealis, FK Comae and Thorne-Zytkow.* `draw_exotic_star`: an old giant (3 Gy or more) is one of them at the catalogue's Galaxy count divided by the old giants of the Galaxy (`OLD_GIANT_SYSTEM_SHARE`, 0.5% of systems, a stand-in). R CrB and FK Comae are flagged merger products; a Thorne-Zytkow object is a single member of kind `thorne-zytkow` flagged `theoretical`, with no companion.
- *Generator epoch.* The period law, the gates and the new draws change generated output; epoch recorded in `tests/golden/epochs.json`.

### 13.3 White-dwarf systems

| Type | Observed (default) | Theory and confidence | Where | Generated and stored | Star model |
|---|---|---|---|---|---|
| White dwarf with a companion | 25% to 26% (32% in the 20 pc sample) [S]; default 30% with 22% wide (`multistar` 7.3) | high | old and intermediate; halo, clusters | F | none |
| Close double white dwarf | **9.5% of white dwarfs** within 4 AU [S, Maoz et al.] | high | old | F: conditional fraction of a white dwarf; axis log-uniform 0.001 to 4 AU; Peters rule; two `white-dwarf` members | none |
| Double white dwarf merger product | **about 10% of white dwarfs** form the high-mass bump; 26% to 54% of wide DWDs have a merger-product heavier star [S] | medium | old | F: 10% of white dwarfs get a mass from the bump (0.8 to 1.1 solar masses), flagged `merger-product` | white-dwarf mass bump |
| White dwarf plus main-sequence star, post-common-envelope pair | most close white-dwarf binaries with a star are post-CE [R]; 22% wide, never interacted | medium | old | F: `interaction` = `post-common-envelope`; axis 0.005 to 2 AU, log-uniform | none |
| Cataclysmic variable (all) | **4.8e-6 pc^-3** from a 150 pc Gaia sample, 77% complete [S, [R3412]]; magnetic types 36% [S]. Per system: 3.4e-5; **1.0e7** in the Galaxy [C]; default 6.5e-4 of white dwarfs | high for the density; population synthesis predicts more [S] | old and intermediate; thick disk | F: an `accreting` pair, orbital period 1.3 to 10 h (donor 0.1 to 0.8 solar masses). Subtypes (polar, intermediate polar, dwarf nova) by weights 36% magnetic | `accreting` state, accretion luminosity (GEN.130 slice c) |
| AM Canum Venaticorum binary | **5.5e-7 pc^-3** (+-3.7) [S, Rodriguez et al.]; 3.9e-6 per system; **1.2e6** in the Galaxy [C]; default 7.5e-5 of white dwarfs | medium | old | F: ultracompact pair, period 5 to 65 min, white dwarf plus helium donor of 0.01 to 0.1 solar masses | kind `stripped-helium` donor |
| Nova, symbiotic, Type Ia progenitor | symbiotic: 276 confirmed [S]; predicted 3e3 to 4e5 by older work, 1.2e3 to 1.5e4 by population synthesis [S]. Default **5e-6 of giants** (4.2e3 in the Galaxy) [I, low]; barium and CH stars about 1% of G and K giants [R] | symbiotic: low; Type Ia double-degenerate merger rate 9.7e-12 per white dwarf per year, 4.5 to 7 times the Type Ia rate [S] | old giants; bulge | F and S (the luminosity-pass giants): a giant with a white dwarf companion 1% (barium) and 5e-6 accreting (symbiotic). A Type Ia explosion is an event and out of scope | `accreting` flag on the white dwarf; barium abundance flag |
| White dwarf plus brown dwarf | rare; part of the brown-dwarf companion rate [R] | low | old | F | brown dwarf kind |
| White dwarf with a planet or debris | 25% to 50% of cool white dwarfs polluted [S]; WD 1856+534 b [S] | medium | old | F (planets, GEN.253) | pollution flag |

**As built (GEN.271).** The white-dwarf binary types are drawn inside `multistar.draw_member_tree` for a placed white dwarf, so the scatter and GEN.631's fill read the same tree. Nothing is stored beyond the member flags and the interaction the tree already had (the new bits `MEMBER_MAGNETIC` 2048, `MEMBER_BARIUM` 4096 and `MEMBER_SYMBIOTIC` 8192 are in `scatter_members.flags`, a SMALLINT).

- *The odds of a companion.* `MEMBER_STATE_BY_PLACED["white-dwarf"]` is now a star at 20.5% and a white dwarf at 9.5% (30% in all, as before; the white dwarf share is the Maoz 9.5% of close double white dwarfs, so every white dwarf companion is one of them). The scatter's bulk screen still draws only the state; `draw_white_dwarf_binary` then names the type of a star companion from one more uniform, at the rows' rates per white dwarf over the 20.5%: AM CVn 7.5e-5 (old population only), cataclysmic variable 6.5e-4 (not young), post-common-envelope pair `post_ce_fraction` times the white dwarf companions (10% of 30%; not young); the rest stay on the period law.
- *Orbits.* A leaf carries an `orbit_au` span; a pair of two single members, one with it, has its axis log-uniform over the span in place of the period law (the stability screens, the Roche limit and the retries are unchanged). Close double white dwarf: 0.001 to 4 AU, and the Peters rule (`survives_gravitational_waves`, now applied to every pair of two white dwarfs) removes the pairs that would have merged by the age. Post-common-envelope: 0.005 to 2 AU; 5% (`HOT_SUBDWARF_POST_CE_SHARE`, a stand-in) of those have a stripped helium star (`stripped-helium`, 0.3 to 0.5 solar masses) in place of a main-sequence star (the table's hot subdwarf row).
- *Cataclysmic variables and AM CVn.* The donor fills its lobe by `ACCRETING_FILL_FACTOR` (1.05), so the pair classifies `accreting` by the GEN.269 geometry: a main-sequence donor of 0.08 to 0.7 solar masses (a stand-in for the donor range; a draw whose period is outside 1.3 to 10 h is made again) and the white dwarf's `magnetic` flag at 36%; an AM CVn donor is a degenerate helium star of 0.009 to 0.13 solar masses (log-uniform; radius 0.03 solar radii at 0.05 solar masses, scaling as the mass to the power -1/3, a stand-in) with periods of 5 to 65 minutes. A white dwarf pair is `accreting` only through these two draws (the builder redraws a white dwarf pair that the period law or a span would put inside its lobe), so the types hold their rates: 6.5e-4 and 7.5e-5 of white dwarfs, measured over 300,000 seeded trees (195 and 23 pairs).
- *Merger products.* `draw_white_dwarf_merger`: 10% of old white dwarfs (knob `merger_product_fraction`) carry `MEMBER_MERGER_PRODUCT` and a rejuvenated `effective_age_yr`. The 0.8 to 1.1 solar mass bump is not drawn here: a white dwarf's mass comes from the star model before the tree, and the tree keeps member 0 at that mass (open; filed).
- *Barium and symbiotic giants.* `draw_giant_companion`: a giant born under 3 solar masses (`BARIUM_GIANT_MAX_INITIAL_MASS_SOL`, a stand-in for G and K) is symbiotic at 5e-6 and barium or CH at 1%; it gets a white dwarf born from a heavier star (`_remnant` with the giant's mass as the lowest progenitor) that has had time to form, on a log-uniform period of 80 to 10,000 days (barium) or 200 to 2,000 days (symbiotic), and carries the flag. The seventh of the star-type uniforms (`STAR_TYPE_DRAWS`) is taken for this.
- *Generator epoch.* The new draws and the changed odds move generated output; epoch recorded in `tests/golden/epochs.json`.

### 13.4 Compact-object systems

| Type | Observed (default) | Theory and confidence | Where | Generated and stored | Star model |
|---|---|---|---|---|---|
| High-mass X-ray binary (supergiant or Be) | 165 known (a 2022 catalogue) [S]; 64 Be/X-ray binaries in 2009, 42 confirmed neutron stars [S]; scale height 150 pc [S]. Default: the 0.3% compact companion of a living massive star (4.2), 40% of which accrete (Be/X-ray 60% of those) | medium | young disk, arms | S and F; `accreting` pair; period days to years | `accreting`, X-ray luminosity |
| Low-mass X-ray binary | 349 known [S]; peaks at the bulge [S]; black-hole low-mass binaries predicted at 1e2 to 1e5 and 1.2e5 black hole plus main-sequence binaries [S]. Default: neutron-star or black-hole companion of a low-mass star 5e-5 of neutron stars and 1e-3 of black holes (then an `accreting` fraction of 1%) [I, low] | low | bulge, old disk | F; period hours to a day | `accreting` |
| Microquasar and ultraluminous X-ray source | microquasars: several dozen [R]; ultraluminous sources are a super-Eddington accretor, mostly extragalactic [R] | low | young | F, a flag on an `accreting` pair; jets along the spin axis (GEN.255) | jet flag |
| Pulsar with a white dwarf (millisecond pulsar) | 78% of millisecond pulsars have a companion (`multistar` 7.3); black widows and redbacks: 41 field spiders, about 80% of binary millisecond pulsars under 2 ms [S]; the companion mass gap 0.07 to 0.1 solar masses is empty [S] | high | old, bulge, clusters | F; companion by class: white dwarf 80%, black widow (0.01 to 0.05 solar masses) 8%, redback (0.2 to 0.4) 4%, neutron star 5%, other 3% (replaces the earlier 10% "other") | `stripped-helium` donors; spiders carry an ablation flag |
| Pulsar with a main-sequence or Be star | PSR B1259-63 (B2Ve, 3.4 yr, e high) and three others with massive companions [S] | medium | young | F; part of the young neutron star's 1% to 3% companion | none |
| Double neutron star | **eight merging in a Hubble time, about twenty known** [S, R]; stock of 6.8e4 to 4.7e5 in the Galaxy [C]; default 0.04% of neutron stars | high that they exist; rate uncertain | all | **S** (promoted pair) and F | none |
| Neutron star with a black hole | none known in the Galaxy; gravitational-wave rate 9 to 86 per Gpc^3 per year [S]; default 0.05% of black holes | medium | all | S and F | none |
| Black hole with a black hole | none known in the Galaxy; 13 to 26 per Gpc^3 per year [S]; default 0.2% of black holes | medium | all; more at low metallicity (halo, clusters, the Clouds) | S and F | none |
| Dormant black hole with a star | **Gaia BH1, BH2, BH3** [S]; 0.1% of black holes with a luminous companion at AU scales [S]; default 1.0% of black holes | medium | disk and halo | S and F | none |
| Dormant neutron star with a star | **Gaia NS1** and 21 astrometric candidates near 1.4 solar masses [S] | medium | disk, halo | S (rare) and F; part of the young-neutron-star companion share (1% to 3%) | none |
| Compact triples and quadruples | PSR J0337+1715, V404 Cygni [S]; about half of black-hole binaries are in triples in a population model [S]. Default 30% of compact binaries have a third member | medium; a quadruple with two compact pairs: low (probability under 1e-6) | all | S and F, same tree; at most three compact members | none |
| Black hole or neutron star in a cluster | NGC 3201: three detached black hole candidates (180 d, 4.3 to 7.7 solar masses minimum) and about 100 black holes in the cluster [S]; globular clusters hold **340 millisecond pulsars**, about half of all known [S] | high | globular clusters, young clusters | F (clusters are fill-time, GEN.159 to GEN.164); cluster rates replace the field fractions: black hole with a star 10%, millisecond pulsar with a companion 80% with a higher neutron-star share, dynamical black-hole pairs from the cluster density | cluster flag on the system |
| Intermediate-mass black hole with a star | candidates only [R] | low | clusters, satellites | S, 1% (section 4.3) | none |
| Magnetar in a binary | none confirmed [R] | low | young | F, 0 by default | none |
| Supermassive black-hole binary | OJ 287 and others (extragalactic) [S] | high | other galaxies, nucleus of a merger | not a star system; the nucleus text flag of `exotic-environments...` 2.4 | none |

**As built (GEN.272).** The accreting compact binaries are drawn inside `multistar.draw_member_tree`, so the scatter and the fill read the same tree. Nothing is stored beyond the member flags and the node flags the tree already had (the new member bits `MEMBER_MICROQUASAR` 16384 and `MEMBER_ABLATING` 32768 fill the SMALLINT of `scatter_members.flags`; the new node bits `NODE_WIND_FED` 4 and `NODE_DISC_FED` 8 are in the TINYINT of `scatter_nodes.flags`). The X-ray luminosity is derived and never stored (`MemberTree.x_ray_luminosities`, `physics/accretion.py`), as are the Be/X-ray, spider and dormant pairs (`be_xray_nodes`, `spider_nodes`, `dormant_nodes`). GEN.952: the members keep their radius (`scatter_members.radius_au`), so a page and the API give an accreting pair its X-ray luminosity, worked out on read by the same call `MemberTree.x_ray_luminosities` makes; it is still derived, not stored.

- *Accreting is explicit.* A pair of a neutron star or a black hole with a single member is `accreting` only through the draws below (the builder redraws a pair that the period law would put inside the donor's lobe, as it does for a white dwarf), so each type holds its rate. `draw_accreting_compact` is the one place the shares are read.
- *High-mass X-ray binaries.* The compact companion of a living massive star keeps its 0.3% (0.2% neutron star, 0.1% black hole); the tree draws three more uniforms for every star (`STAR_TYPE_DRAWS` is now 10): 40% of those companions accrete (knob `accreting_fraction`), and 60% of the accreting ones are Be/X-ray binaries (the donor flagged `be`, a decretion disc, period log-uniform over 10 to 300 days on the eccentric thermal law, young population only), the rest are wind-fed supergiant binaries (1.4 to 100 days). A period range whose closest axis would put the donor inside `MAX_FILL_FACTOR` of its lobe starts at the closest allowed axis and reaches at least four times past it.
- *Low-mass X-ray binaries.* They are part of the `star` odds of a neutron star and a black hole, not an addition to them, so a black hole's 1.25% total stays: 5e-5 of neutron stars and 1e-3 of black holes over those odds (0.25% and 10% of the star companions, `draw_low_mass_xray_pair`; not a recycled pulsar and not the young population) have a low-mass donor (0.1 to 1.5 solar masses) on the interacted period law, and 1% of those accrete: the main-sequence donor fills its lobe at a period of 4 to 24 hours (`ACCRETING_FILL_FACTOR`), by `overflow`.
- *Microquasars.* An accreting neutron star or black hole of either kind launches jets at `MICROQUASAR_SHARE_OF_XRAY_BINARIES` (10%, a stand-in from the several dozen microquasars against the 165 high-mass and 349 low-mass X-ray binaries): the flag `microquasar`, set on the accretor.
- *Luminosity.* `L = eta Mdot c^2`, at most the accretor's Eddington limit (`4 pi G M m_p c / sigma_T`). A neutron star's `eta` is its surface potential, a black hole's `ACCRETION_EFFICIENCY_BLACK_HOLE` 0.1. A wind-fed pair captures the Bondi-Hoyle share of the donor's wind (loss rate 1e-6 solar masses a year at 1e5 solar luminosities, to the power 1.5 of the luminosity from `1.5 M^3.5`; speed 1.5 times the escape speed), a Be/X-ray pair the same from the disc (1e-8 solar masses a year at 20 km/s) at the closest approach, a lobe-filling donor its mass over 5e8 years (a main-sequence donor) or 1e7 (a giant). Over 40,000 seeded massive stars the median wind-fed pair shines at 4e36 erg/s and none above the limit.
- *Millisecond pulsars.* The recycled pulsar's table now has the rows' classes as states: white dwarf 66.1%, black widow 6.24% (0.01 to 0.05 solar masses), redback 3.12% (0.2 to 0.4), other 2.34% (a star of 0.1 to 1 solar mass) and neutron star 0.2% (78% with a companion in all, as before). Spiders have periods of 0.1 to 1 day, are `detached`, and the pulsar carries `ablating`; the 0.07 to 0.1 solar mass gap stays empty because the two mass ranges do not reach it. The row's 5% neutron-star share is not drawn: the design note's 4.4 budget (the 0.04% of double neutron stars) leaves 0.2% of recycled pulsars, and five percent of them would make about ten times the double neutron stars the gravitational-wave rate allows; the white dwarf share took the difference.
- *Dormant companions (Gaia kind).* A `detached` pair of a neutron star or black hole and a main-sequence star or giant, on the interacted law (0.3 to 1e4 days); `dormant_nodes` lists them.
- *Compact triples.* The `compact-triple` row's 30% now applies to neutron stars as well as black holes (a third member, a compact object of either kind in 10% of them, else a star; at most three compact members). The quadruple with two compact pairs (theoretical, under 1e-6) is not drawn, and magnetar binaries stay at 0.
- *Generator epoch.* The new draws, the pulsar table and the rule that a neutron star or black hole pair never accretes by geometry alone move generated output; epoch recorded in `tests/golden/epochs.json`.

### 13.5 Substellar companions, planetary-mass companions and disks

| Type | Observed (default) | Theory and confidence | Where | Generated and stored | Star model |
|---|---|---|---|---|---|
| Brown dwarf close companion | under 0.6% inside 3 to 10 AU for Sun-like stars (the "desert") [S] | high | everywhere | S and F: kind `brown-dwarf` member; GEN.252 | L, T, Y states |
| Brown dwarf wide companion | **3.2%** at 28 to 1,590 AU (6.8% in one thesis) [S] | medium | everywhere | S and F, default 3.2% | as above |
| Brown dwarf pair or brown dwarf plus low-mass star | about 15% of brown dwarfs are in pairs [R] | medium | everywhere | F | as above |
| Planetary-mass companion (a few Jupiter masses, wide) | directly imaged companions are a few per cent of young stars at tens to hundreds of AU [R] | low | young | F, a member of kind `brown-dwarf` below 13 Jupiter masses, flagged `planetary-mass` | none |
| Circumbinary planet | **at least 10%** of close eclipsing binaries host a planet above 6 Earth radii at periods under 300 d if coplanar (47% minimum if isotropic); giants are rare [S] | medium | everywhere | F only (planets); Holman-Wiegert P-type zone from the node | none |
| Circumbinary disk | protostar binaries, post-AGB binaries, planetary-nebula binaries: disks but no planets [S] | medium | young, evolved | F, a descriptor on the node, no bodies (`exotic-environments...` 1.3) | `disk` |
| Planet in a wide pair (S-type) | planets known in wide binaries [S, `multistar` 5.6] | high | everywhere | F only | none |
| Rogue planet pair (a rogue and its moon, or two rogues) | free-floating pairs reported [R, disputed] | low | clouds | F, existing moon rule | none |

### 13.6 Ejection and runaways

| Type | Observed (default) | Theory and confidence | Where | Generated and stored | Star model |
|---|---|---|---|---|---|
| Runaway star | **10% to 30% of O stars and 5% to 10% of B stars** [S, [R1558]]; a bow-shock survey finds 23% [S] | two-thirds binary-supernova, one-third dynamical for massive runaways [S, Dray et al.] (disputed) | young disk, clusters | the existing runaway flag; new `runaway_origin` (supernova or dynamical) and, for the supernova channel, a **compact companion left behind** at the rate the survival rule gives | `runaway_origin` |
| Hypervelocity star | 1e3 to 1e4 [S] | high | halo | existing scatter; single (the Hills mechanism) | none |
| Companion unbound by a collapse | 83% of non-interacting massive binaries are unbound by the first collapse [C] | high | young | not stored; the survivor is an ordinary row, with a flag if it is a runaway | none |

**As built (GEN.327).** `multistar.draw_young_group` draws the trapezium of the catalogue (`young-group`, 5% of the O and early-B stars, birth mass 9 solar masses and up, younger than `YOUNG_GROUP_AGE_MAX_YR` = 0.3 Myr), and `draw_member_tree(..., groups=True)` calls it first; only the bright-star scatter asks (`scatter_companions.star_tree`), because a system build holds pairs, not a group. A star that is not eligible draws nothing, so its tree is as before. The group has 3 to 5 stars and a `decay_age_yr` drawn log-uniform over 1e4 to 4e4 years. Younger than that, the tree is the group: one `group` node (separation = the group radius, period = the crossing time at the dispersion, no eccentricity) and its members, each with an `offset_au` from the centre and a `velocity_kms` (each component a draw of the dispersion), shifted so the centre of mass is at rest at the origin, and the group's decay age (the sixteen bits of the flags column are all taken, so a group member has no flag of its own: the node's kind and the decay age say it has not decayed). Older than that, the tree is the result: the lightest star other than the placed one is ejected, and the rest are drawn as a hierarchy no wider than `YOUNG_GROUP_REMNANT_MAX_AXIS_AU`. The ejected star is on `MemberTree.ejected`, flagged `runaway_dynamical`, and is **not** a row of its own: the runaway rows the scatter draws (`draw_runaway_origin`, 20% of O and 7.5% of B stars, a third by the dynamical channel) already carry the catalogue's share of ejected stars, and a row for each ejection would count them twice. Since 4e4 years is short against 0.3 Myr, most of the 5% leave a remnant and a stored group is rare. `scatter_members` stores the offset, velocity and decay columns (they were in the schema since migration 0104; the writer now fills them); a stored group is not built by a system build (`stored_trees.plan` takes only pair nodes, so the placed star is built alone as for any tree it cannot hold) until GEN.264 builds the stored members. The uncharted star page lists the group with the life it has left (decay age less the star's age), and `scatter_systems.groups_past_decay` is audit 14. The member count weights, the radius range, the dispersion (2.5 km/s) and the remnant axis cap (300 AU) are stand-ins (kind P) in `tuning.py`, each with its bracket: replace them when a count by size, the Trapezium's radius and its measured dispersion are read. GEN.328 draws the same group at fill (below).

**As built (GEN.328).** The fill draws the same group (`nebula_contents.member_group`, on a stream of the member's own seed, so every sector that reaches the member draws the same one) for each member of an emission region (D or E) or a field cloud that is an O or early-B star younger than `YOUNG_GROUP_AGE_MAX_YR`: `multistar.draw_young_group` at the `young-group` share. The protostars of a core are all under `PRE_MAIN_SEQUENCE_MAX_MASS_SOL` (8 solar masses), so a core holds no group, and the ages of an emission region (a uniform draw up to 10 or 20 Myr) leave only a few in a hundred members under 0.3 Myr: groups are rare at fill by the numbers of the catalogue, and a test forces the age to measure the share. A system build cannot hold a `group` node (`stored_trees` pairs, GEN.264 for the stored members), so a group that has not decayed is built as separate systems: the member keeps the point it was given, each mate is an ordinary young system of its own birth mass and the member's age at the point plus the difference of the two offsets (`run_sector.add_group_stars`), and each sector builds the stars whose points fall in its cube, so a group across a border is made once. The entries carry `group_seed` (the member's seed names the group), `group_decay_age_yr` and `group_offset_au`; they are in memory, as `system_seed` is, and the velocities of the group's tree are not used (a fill system holds no velocity of its own). A group older than its decay age is the hierarchy it left (`stored_trees.drawn` of the tree without its brown dwarfs and remnants) and marks no group; its ejected star is not a row, as at scatter. An emission region's reach (`emission_reach_pc`) grows by two radii of the widest group (`group_reach_pc`, 0.29 pc) so a neighbouring sector still finds a mate.

**As built (GEN.329).** The catalogue's rates are read from `SYSTEM_TYPE_RATES["runaway-star"]` and `["unbound-by-collapse"]`, and `multistar.draw_runaway_origin(initial_mass_sol, runaway_unit, origin_unit, bound_unit)` is the generator they name. An O star (birth mass 16 solar masses or more) is a runaway with the default 20% and a B star (2.1 to 16) with 7.5%, inside the catalogue's 10% to 30% and 5% to 10%; any other star has no observed rate and keeps the field roll of the fill. A runaway's `origin` is `supernova` for two thirds of them and `dynamical` for the rest, and a supernova runaway keeps a compact companion bound with the share the collapse leaves bound, 1 - 0.83 = 17%; the unbound 83% are not stored, the survivor is an ordinary row. The three uniforms are derived from the row's seed (`scatter_companions.runaway_draw`, three more salts of `screen_unit`), so the scatter that stores the star's tree and the fill that builds the star agree without storing the answer for a single star. A runaway's tree carries the flag of its channel on member 0 (`MEMBER_RUNAWAY_SUPERNOVA` 128, `MEMBER_RUNAWAY_DYNAMICAL` 256, in `scatter_members.flags`), and a supernova runaway that kept its companion has a neutron star or a black hole in it (`draw_member_tree(bound_runaway=True)`: the table of its class decides which, the progenitor is heavier than the runaway and has had time to collapse by the system's age; a runaway too young for any progenitor, such as a 20 solar mass star at 5 Myr, finds one in about a quarter of the draws). The fill (`run_sector._flag_scattered_runaway`) makes the system of a pre-placed O or B star a `runaway` from the same answer, its speed and direction drawn from a stream of the row's seed, and marks it decided so the field roll (`flag_fast_stars`) leaves it alone; the field roll gives the runaways it makes an origin too, by the same split. `StarSystem.runaway_origin` (serialization version 5) is stored in the flags of the runaway's first `stars` row (`insert_star(flags=...)`, the bit of its origin) and read back by `load_star_system`; a hypervelocity star has none and stays single (no tree is drawn for it). No schema change: `stars.flags` is the v104 column. Generator epoch 56.

### 13.7 Tuning the emergent types

Contact, Algol, post-common-envelope and accreting fractions come out of the period law, the Roche test and the evolution, and the first build will not match the observed fractions. The tuning dict has one knob per family, all against the audit counts in section 15:

- **tight-tail weight** (default 6.0 as built, GEN.270; 1.0 is the plain draw): scales the share of pairs with a period under 1 day of primaries under the massive limit until the contact share is 0.2% of stars and the Algol share 0.1% to 0.2%;
- **post-CE fraction of white-dwarf companions** (default 0.10, the Maoz 9.5% for double white dwarfs and a matching value for white dwarf plus star);
- **accreting fraction of compact companions** (default 1% to 40% by class, rows above);
- **merger-product fractions** (massive 8%, field 0.5%, white dwarf 10%).

## 14. What changes in the star model and the stored form

**Star model** (`physics/stellar_evolution.py`, `Star`, `star_labels.py`; these are also the TODO items 13 to 16 below):

1. **New `kind`s**: `pre-main-sequence` (GEN.242), `protostar`, `brown-dwarf` with L, T, Y states (GEN.252, built: `brown-dwarf-state.md`), `stripped-helium` (hot subdwarf, helium star, helium donor), `thorne-zytkow`, `wolf-rayet` (GEN.243), `white-dwarf` subtypes (normal, polluted, merger-product, helium) and the compact kinds already there.
2. **Interaction state on a pair** (`scatter_nodes.interaction` and `orbit_nodes.interaction`): `detached`, `semi-detached`, `contact`, `accreting`, `post-common-envelope`, `merger-product`.
3. **Per-member flags** (`scatter_members.flags`): `mass_reversed`, `mass_gainer`, `merger_product`, `be`, `runaway_origin` (none, supernova, dynamical), `fragile` (ultra-wide), `theoretical`.
4. **Roche geometry and accretion luminosity**: Roche-lobe radius (Eggleton formula) and filling factor, mass-transfer luminosity `L = eta Mdot c^2` for accreting pairs (GEN.130 slice c), Eddington limit for the accretor.
5. **Rejuvenated age** for a merger product (an effective age field; the real system age stays).
6. **Group nodes for non-hierarchical systems**: a node of kind `group` with N members, each with an offset and a velocity at an epoch; its internal motion is a straight line until the decay time, then the generator's remnant. The nightly update does not integrate it (the system is young; at the decay time it is replaced by the remnant at the next fill or edit).

**Stored form changes** to section 5: `scatter_members.slot` and `scatter_nodes.slot` become SMALLINT (a group has N slots); `scatter_nodes.kind` gains `group` and `compact`; new nullable columns `interaction` (VARCHAR(20)) on nodes, `state`, `flags` (SMALLINT) and, for group members only, `offset_x_au`, `offset_y_au`, `offset_z_au`, `velocity_x_kms`, `velocity_y_kms`, `velocity_z_kms`, `decay_age_yr` on members. One column per component; no JSON. A group's lifetime is `decay_age_yr` and the audit rejects a stored group older than it.

**The dict.** `SYSTEM_TYPE_RATES` has one entry per row of 13.1 to 13.6: the reference class, the default, the low and high of the range, the evidence tag, the populations that may host it (the profile's `mix`), and the generator function. A completeness test (the same one the registry has, GEN.250) checks that every row of this catalogue has an entry, a generator and an audit count.

## 15. Audit extension for the catalogue

Extends section 9:

10. **Type counts.** Over a seeded fill of at least 1e6 systems per region class, each type's count lies within four Poisson sigma of the default (or, where the default has a range, within the range for the sample size): contact 0.2% of FGK, Algol 0.1% to 0.2%, cataclysmic variable 3.4e-5 of systems, AM CVn 3.9e-6, double white dwarf 9.5% of white dwarfs, brown dwarf wide companion 3.2%, runaway 10% to 30% of O stars.
11. **Region gates.** No type appears where its gate forbids it (a Be/X-ray binary in the halo, a trapezium older than its decay time).
12. **N tail.** The number of systems with N = 8 to 12 matches the tail table (1.1e6 to 1.1e2 in the default galaxy) within four sigma on a sample scaled to it.
13. **Interaction consistency.** Every `semi-detached` pair has a donor whose radius is at least its Roche lobe at the system age; every `contact` pair has both; no `accreting` pair has a non-compact receiver unless it is a mass gainer.
14. **Group decay.** Every stored `group` is younger than its `decay_age_yr`; a drawn group older than it is replaced by its remnant.

## 16. Decisions that this section changes

- **The cap is 12, not 7.** Evidence: ten candidate systems with multiplicity up to 12 [S, [R1529]]; theory several hundred [S]; stability allows 12 or more levels (3.3); a septuple protostellar system found in 2025 [S]. At the tail table's ratio the Galaxy holds 1.1e6 eight-star systems and 1e2 twelve-star systems [C]. The earlier note's 7 was the **confirmed** ceiling; 12 is the highest *proposed*. The prevalence knob (GEN.129) scales the N >= 5 columns, and a directive forces N.
- **Short-lived young groups are in.** The earlier "trees only" default (decision 4) becomes: stored systems are hierarchical, **except** young non-hierarchical groups under their decay age (13.1), drawn at fill in star-forming cores and, rarely, at scatter. A group older than its decay age is never drawn; its remnant is.
- **Types are mostly emergent.** The classifier of 13 produces contact, Algol, post-common-envelope and accreting pairs from geometry and evolution; the tuning dict aligns them to the observed fractions.
- **Observed rates are defaults; theoretical rates are marked.** A type with no observation (Thorne-Zytkow, magnetar binary, quadruple with two compact pairs, planetary-mass companion beyond a few per cent) takes the predicted rate, carries the flag `theoretical`, and the page says so.

## Evidence notes

[R] items to confirm when paper access is allowed: the Kelvin-Helmholtz time; the number of known high-mass X-ray binaries (about 130); the luminosity and mass of Cygnus X-1; the double pulsar's period and merger time; the abundance of barium stars among giants (about 1%); the planetary-nebula central star close-binary fraction (20%); that Hills ejection leaves a hypervelocity star single; the density of Milky Way-equivalent galaxies used to turn gravitational-wave rates into per-galaxy rates (0.0116 per Mpc^3, 1.16e7 per Gpc^3); the [R2539] count of 1.2e6 black-hole pairs (read in a search-result summary of [R3348], not at source); the Gaia 0.1% prediction (a summary of [R3233]); every Gaia and pulsar orbital figure.

Added by the catalogue: the symbiotic rate (5e-6 of giants, low), the field merger-product rate (0.5%), the massive over-contact rate (1%), the trapezium share (5% of young massive systems) and the Thorne-Zytkow count (100) are inferred, not observed; every Roche-classifier fraction depends on the tuning of 13.7.

Weakest numbers: the black-hole pair, black-hole with neutron star and black-hole with star fractions (a range of 80, 90 and 20 times); the direct-collapse mass that drives the wide-channel Monte Carlo; the 30% tertiary share; the compact-companion share of massive stars (0.3%). All are single entries in the dict, and a changed value changes the next scatter only.

[C] values come from scripts kept in the session scratchpad and not shipped: `counts.py`, `rows.py` (row counts and class mix from the repository's own `star_population` and `phenomenon_scatter` functions), `mc.py` and `mc2.py` (the survival Monte Carlo: Sana period law on 0.15 to 5.5 in log days, 90% of massive primaries with a companion, mass ratio uniform 0.05 to 1, thermal eccentricity above 7 days, instantaneous collapse on a random orbital phase with an isotropic Maxwellian kick, a system counted as interacting and excluded when its pericentre is under 9.3 AU, twice the 1,000 solar-radius giant phase), and `stock.py` (Peters axes and the stock of compact pairs from a merger rate: with a delay-time distribution proportional to 1/t between 50 Myr and 13.8 Gyr the unmerged tight pairs number the merger rate times 13.8 Gyr / ln(13.8 Gyr / 50 Myr) = 2.5 Gyr, and the share of a log-flat 0.003 to 100 AU draw that merges in a Hubble time is 0.19, 0.26 and 0.33). The Monte Carlo is a size check, not a population synthesis; it has no mass transfer, common envelope or tides.

## Sources

Seen in search-result summaries (no paper was opened): [R0006, R1076, R1105, R1239, R1381, R1709, R2407, R2463, R3233, R3557, R4039, R4441]; [binary black hole formation and metallicity][R4512]; [R3108, R3333]; [triple evolution (Toonen et al.)][R4404]; [R0074, R0075, R1363, R4203, R4279, R4386]. Repository files read: `generation/phenomenon_scatter.py`, `bright_stars.py`, `star_population.py`, `run_sector.py`, `system.py`, `binary.py`, `wide_binary.py`, `tuning.py`; `docs/reference/database-schema.md`, `docs/TODO.md` (GEN.129, GEN.130, GEN.223, GEN.250, GEN.251); the design notes linked above and `phenomenon-scatter-mass-cut.md`.

Added for the catalogue (search-result summaries; no paper opened): cataclysmic variable space density from a 150 pc Gaia sample [R3412]; AM Canum Venaticorum space density (Rodriguez et al.); close double white dwarfs, 9.5% of white dwarfs (Maoz et al.); contact binary frequency and tertiaries (Rucinski; Tokovinin et al.; Kepler eclipsing-binary morphology); blue stragglers from triple mergers (Perets and Kratter; Naoz and Fabrycky); symbiotic star counts and predictions (Merc et al.; Lu et al.); Thorne-Zytkow candidates HV 2112 and HV 11417 (Beasor et al.); X-ray binary catalogues (Liu et al.; Avakyan et al.); spider millisecond pulsars; [R1529] on high-order multiples; the Nu Scorpii septuple and the NGC 6334IN protostellar septuple (2025); trapezium-type multiples (Allen et al.; Orion Trapezium decay); runaway O and B stars ([R1558]; Dray et al.); circumbinary planet occurrence (Armstrong et al.; Martin and Fabrycky); brown-dwarf companion statistics. Where a figure here has no link, it is a snippet read in a search, tagged [S], and the weakest ones are listed under Evidence notes.
