# Multi-star systems at scatter time: census, limits, probabilities and the stored tree

Boss (2026-10-15, 5:51 PM EDT): probabilities accurate but "no forced representation in a sector (except if one landed there in a scatter before generation)"; stars and star-like objects (former stars, protostars) are generated, and "if a star, black hole, neutron star, etc has a companion I want it generated at scatter time", including neutron stars with neutron stars, black holes with neutron stars, black holes with black holes and stars with neutron stars or black holes, because "one of the points of our scattering is to represent gravity and such systems would definitely impact gravity". He asked for research into how many multi-star (or star-like object) stable systems exist, the theoretical limits, and the probability of each, then an expansion of our multi-star definitions.

This note does that. It extends [multistar-and-compact-systems.md](multistar-and-compact-systems.md) (the stability criteria, the 2 to 7 star hierarchy, the compact-primary rules; GEN.128 to GEN.130), [phenomena-catalogue-and-fill.md](phenomena-catalogue-and-fill.md) (the kind list and the fill stages; GEN.250, GEN.251) and [exotic-environments-planets-and-compact-binaries.md](exotic-environments-planets-and-compact-binaries.md) (Peters time, V404 Cygni, survival of wide companions). It does not repeat their numbers; it says where this note changes them (section 10).

Informs: GEN.129, GEN.130, GEN.223, GEN.109, GEN.250, GEN.251, GEN.100, GEN.242, GEN.243, GEN.252, and the new items in section 11.

Status: design, 2026-10-15; nothing here is built. Where Boss's words leave a gap the note takes a default, marked "Default" and listed in section 12.

Evidence tags: **[S]** seen in this repository or in a search-result summary (the papers themselves could not be opened: arXiv is unreachable from this environment, so every literature number is a snippet, never a table read from the paper); **[C]** computed in this session (scripts named in "Evidence notes"); **[R]** recalled, unconfirmed; **[I]** inferred, with the step shown. A figure keeps the tag of its weakest input. Times shown to Boss are US Eastern (EDT) only.

## Summary

- **The rule changes in one place.** Boss's 5:39 PM EDT rule said scatters carry no systems, planets or moons. His 5:51 PM EDT message narrows that: a scatter row still carries **no planets, moons or belts**, but it **does carry its star-like companions** (stars, brown dwarfs, white dwarfs, neutron stars, black holes, pre-main-sequence stars), generated when the scatter is drawn. Planets, moons and belts stay fill-only.
- **A scatter row becomes a system key.** Its point is the system's barycentre, its mass is the system's total mass, and its members and their orbits sit in two new child tables with one column per component and no JSON (section 5). A fill builds the system from those rows instead of drawing the companions again, so a companion that landed in a scatter is the companion the sector gets.
- **The cost is small.** At the default galaxy the scatter holds about 1.13e6 star rows and 5.9e7 black-hole rows [C]. Companions add about 4.9e6 member rows (the placed object included) and 3.0e6 orbit nodes, roughly 0.6 GB against the 8 GB the base rows already need, and about half an hour of CPU at scatter time (section 6).
- **The census says what the probabilities are** (section 4). Stars follow the existing O to M table: a scattered O or B star has a mean of 3.1 or 2.3 stars (the rows' mix gives 2.77 and 2.28) [S, C]. A black hole has a companion at about 1.25% (star 1.0%, black hole 0.2%, neutron star 0.05%), not the 10% to 30% the earlier note took from three Gaia systems [R, I]. A neutron star has a neutron-star companion at about 0.04% [C, I]. Every default is one number in one tuning dict, with a range.
- **Neutron-star pairs would be missed by the mass cut, so promote the pairs.** A neutron star never passes the 14 solar mass cut, so a double neutron star would be a fill-time object. A pair is rare (about 1.7e5 of them in the default galaxy) [C, I], so the scatter places every bound compact pair (neutron star with neutron star, with black hole, black hole with black hole) whatever its mass, about 3e5 extra rows. Default on; Boss decides (section 4.5).
- **The theoretical limit is high and the observed one is 7.** A hierarchy can in principle hold a dozen levels between a contact binary and the galaxy's tidal limit [C]; seven stars is the most ever proposed with none certain [S]. The default cap stays 7; a directive raises it. Non-hierarchical systems (trapezia) decay in 1e4 to 1e5 years [S], so they are not stored as stable. Compact members add three filters: survival of each collapse, the gravitational-wave merger time, and the relativistic quench of Kozai-Lidov cycles (section 3.3).
- **Motion needs almost nothing new** (section 7). The barycentre point moves with the system's total mass (GEN.223); a member's position is a function of time computed on read from its stored phase, epoch and period, so the nightly update writes nothing per member.
- **Eleven proposed TODO items** (section 11), and five rewordings of existing items. Six decisions for Boss (section 12).

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

- **Rows.** `bright_stars` holds a star (full star record, a `seed`, a point in mpc); `phenomenon_scatter` holds `kind`, `subtype`, a point, a `seed`, `mass_solar`, `map_luminosity_sol` (`phenomenon_scatter.py:73`, `docs/database-schema.md` 1795 and 1948). A row is a point and a seed. A fill later builds the object from the seed.
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

The star numbers are in `multistar-and-compact-systems.md` 3.1 to 3.3 and stand (Raghavan 2010: 56% single, 33% double, 8% triple, 3% higher among FGK dwarfs within 25 pc [S]; Tokovinin 2014 [S]; Moe and Di Stefano 2017: 0.50 companions per solar-type primary rising to 2.1 for O stars [S]; the observational ceiling is seven stars, none certain [S]). New in this note, from search-result summaries [S]:

| Quantity | Value | Source |
|---|---|---|
| Triple fraction of early B and O stars | about 50% and 70% | Moe and Di Stefano 2017, as quoted in Stegmann, Antonini and Moe 2022 |
| Only 16 +- 8% of the stars that dominate the core-collapse supernova rate are single | | Moe and Di Stefano 2017 |
| White dwarfs: share with a companion | 25% to 26% (232 white dwarfs within 25 pc: 74% single, 26% not); the 20 pc sample gives 32% | Holberg et al.; Toonen et al. 2017 |
| Black-widow and redback spiders | 41 confirmed and suspected in the field; with the redbacks they are about 80% of binary field millisecond pulsars with spin period under 2 ms; companion masses 0.07 to 0.1 solar masses are almost absent | arXiv 2210.16295 |
| Surviving the first supernova | binary survival about 14% (+20, -10), neutron-star and black-hole remnants in similar numbers; 1% to 10% of core-collapse supernovae occur with a compact companion present | Renzo et al. 2019 as quoted in arXiv 2401.02389 |
| Fraction of massive-star binaries that end as a merging black-hole pair | about one in eight interacting binaries massive enough to make black holes (COMPAS, low metallicity) | arXiv 2411.02484 |
| Gaia dormant neutron-star candidates | 21 astrometric binaries with solar-type stars and dark companions near 1.4 solar masses; Gaia NS1 (1.90 +- 0.04 solar masses, 731 d, a halo G star) | El-Badry et al. (arXiv 2402.06722, 2405.00089) |
| Gaia black holes with luminous companions at AU scales | about 0.1% of Milky Way black holes (prediction of about 30 in Gaia DR4); up to about 900 detectable in the thin disc | arXiv 2502.03527; Kotko, Banerjee and Belczynski 2024 |
| Black-hole pairs in the Galaxy | about 1.2e6 (cosmological simulation, mean system mass 28 solar masses) | Lamberts et al. 2018 as quoted in arXiv 1908.08775 |
| Double neutron stars: merger rate per Milky Way | 42 (+30, -14) per Myr (Pol 2019, eight merging systems), 37 (+24, -11) (Pol 2020), 32 (+19, -9) (Desvignes 2019); gravitational-wave based 5.3 (+17.6, -4.9) and 11.0 (+19.2, -8.2) per Myr | arXiv 1811.04086, 2107.13307, 2402.04658, 2604.05059 |
| Local merger-rate density (GWTC-4.0, 90% intervals) | binary neutron star 7.6 to 250 (published 8.8 to 250), neutron star with black hole 9.1 to 84 (9.3 to 86), binary black hole 14 to 26 (13 to 26), per Gpc^3 per year | arXiv 2508.18083 |
| Trapezium-type systems decay in | 1e4 to 4e4 years (bound ones); after two crossing times 41 of 100 still trapezia; the ejected star is usually the lightest | arXiv 1809.03537, 1701.03440, 1908.11437 |
| Black-hole binaries in triples | about half of those formed in the TSE model are in triples | Stegmann, Antonini and Moe 2022 |

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

1. **Hierarchy depth.** A level needs an outer-to-inner axis ratio of about 3 or more for equal masses (Mardling and Aarseth 2.9 to 3.3 coplanar [S, `multistar...` 5.1]) and observed period ratios are mostly above 10, i.e. axis ratios above about 4.6. The tightest stellar orbit is a contact binary at 0.01 AU and the widest bound orbit the generator keeps is 1e4 AU (20,000 AU for the root), a span of 1e6 to 2e6 in axis. A chain of axis ratio 3 fits `ln(1e6)/ln(3)` = 12.6 levels, and a ratio of 5 fits 8.6 [C]. So the stability criteria allow a chain of **a dozen or more stars**, and the earlier note's "theory allows several hundred" [S, Gebrehiwot et al. 2016] needs branching. Nobody has seen more than seven. **Default: cap at seven** (`MAX_STARS`), a directive raises it.
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

A member is drawn with an initial mass, then its state follows from the system's age with the same evolution code a single star uses (`physics/stellar_evolution.evolve_star`). A member lighter than the placed object (the usual case: the placed object is the heaviest) lives longer than it does, so it is alive: main sequence, giant, or **pre-main-sequence** if its age is under its Kelvin-Helmholtz time (`t_KH` of about 31 Myr x (1 solar mass / m)^2.5 [R]; this is what makes the companions of a 5 Myr-old O star young, GEN.242), or a **protostar** if the system is under 0.5 Myr (only inside the cores of GEN.241; a scatter row is never that young). A **brown dwarf** companion has the desert rule (under 1% inside 5 AU, a few per cent wide; GEN.252). A **former star** companion appears only through the explicit compact-companion probabilities below, because a companion that is heavier than the placed object at birth contradicts "the placed object is the heaviest"; the compact probabilities stand for the mass-transfer reversal and the heavier-companion cases together. A remnant is allowed only if its progenitor's lifetime is shorter than the system age and the placed object's lifetime is longer (for a living 14 to 25 solar mass star the window is its last third); the table value is the mean over all ages and the generator divides it by the window's share and applies it inside the window.

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
| `slot` | CHAR(1) | `a` or `b`. |
| `kind` | VARCHAR(20), CHECK | `main-sequence`, `giant`, `white-dwarf`, `neutron-star`, `black-hole`, `brown-dwarf`, `pre-main-sequence`, `protostar`. |
| `mass_solar` | DOUBLE | Current mass. |
| `initial_mass_solar` | DOUBLE | Birth mass (the progenitor's for a remnant); a star's current state is rebuilt from it and the age. |
| `seed` | BIGINT UNSIGNED | The member's own build seed (details, planets later). |

**`scatter_nodes`** (N - 1 rows for N members), the same names and units as GEN.129's `orbit_nodes` so a fill copies column for column:

| Column | Type | Notes |
|---|---|---|
| `source`, `source_id`, `node_no` | | Primary key; `node_no` 1 is the root. |
| `parent_node_no`, `slot` | | NULL parent for the root. Children are found by their parent and slot. |
| `mass_solar` | DOUBLE | The subtree's mass. |
| `separation_km`, `eccentricity`, `inclination_deg`, `ascending_node_deg`, `arg_periapsis_deg` | DOUBLE | The pair's orbit (the relative orbit of child b about child a). |
| `mean_anomaly_deg`, `epoch_unix`, `period_years` | DOUBLE | The phase at an epoch. |
| `kind` | VARCHAR(16) | `close`, `wide` or `compact` (for naming and the checks). |
| `flags` | TINYINT | bit 0 interacted (drawn from the tight-channel law), bit 1 `kl_active`. |

One column per component, no JSON. There is no foreign key (a row table of 6e7 entries with a polymorphic parent), so the build deletes the child rows of a built row in the same transaction, and a nightly check (extends `check-db`) counts child rows without a parent.

`phenomenon_scatter_classes` gains no column: a compact pair is counted under its placed object's class, with a `members` count in the Phenomena table.

## 6. The generator and its cost

### 6.1 One function, two call sites

`draw_member_tree(placed, age, rng)` returns the members and nodes. It is GEN.129's hierarchy generator (design 5.1 to 5.3 of the earlier note: draw N from the O to M table, masses, shape, root separation, then each pair top-down under the stability criteria and the Kozai screen) with three additions: the **member kinds** and the compact-companion probabilities of section 4, the **compact filters** of 3.3 (Peters, quenched Kozai, wide-companion rule), and the **state of each member** at the system's age (4.1). The scatter pass calls it for every placed star and black hole, and the fill calls it for every object the fill creates (stage D's other systems, stage G's below-cut remnants), so the two places draw from one set of laws.

### 6.2 Streams and determinism

The tree is drawn from `draw.bound(row.seed)` with the suffix `members`, in a fixed order (GEN.55, GEN.56): the existing rule that a new stage never moves an older stage's draws. The planets later come from the suffix `planets`. A re-plan with the same galaxy seed regenerates the same members; after a plan the stored members are authoritative, so a later change to the laws changes the next scatter, not the existing rows (no drift between a stored summary and a regeneration).

### 6.3 Bulk screening

A black-hole scatter row is made in bulk by numpy (`_band_rows`, about 1 microsecond a row). The companion decision must not turn that into a Python loop over 6e7 rows. Therefore the bulk pass draws **one uniform per row** from the shell's stream and compares it with the cumulative companion probability of the row's class (1.25%); only the hits (about 7.4e5 black-hole rows) call `draw_member_tree`.

### 6.4 Density bookkeeping: a row is a system

The rates in `PHENOMENON_DENSITY_PC3` count **objects**. A black-hole pair would be two objects and, if both were placed independently, a double count, and a pair of which one member is under the cut would be counted once by the scatter and once by the fill. Rules:

- A scatter row is a **system keyed by its placed object**: its rate is the object rate divided by `1 + p_placed`, where `p_placed` is the expected number of **same-class** companions per row (0.20% for a black hole; 0.04% for a neutron star), so the black-hole density is unchanged. The correction is 0.2%, below the rate uncertainty, but it is a rule so the audit can check it.
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

The TODO thread numbers them and assigns the lanes; next free GEN.260, TEST.162, DOC.50, UX.114, MAP.208. Order is build order.

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
11. **Cap and directive** (Core, small): `MAX_STARS` setting (default 7) and a `--stars N` directive, so a user may force N or raise the cap; a prevalence knob for the N >= 5 columns, as GEN.129 already has.

## 12. Decisions for Boss, with the defaults taken

1. **Black-hole companion fraction 1.25%, not 10% to 30%** (4.3). The old number was an inference from three Gaia systems; the censuses support 1.25% (range 0.2% to 7%). Default: 1.25%. The tuning dict allows any value.
2. **Seven stars stays the cap** (3.3). Theory allows a dozen levels. Default: cap 7, directive to raise it.
3. **Promote every bound compact pair to the scatter** (4.5): about 3e5 rows, makes neutron-star pairs exist at scatter time. If off, neutron-star pairs appear only when their sector fills. Default: on.
4. **Hierarchical trees only; no stored trapezia** (3.3). Default: yes; a decaying trapezium in a protostar core is a later flavour item.
5. **Delete a row's members and nodes when its sector builds it** (8) rather than keeping two copies. Default: delete.
6. **Pair rates from the middle of the evidence** (4.3, 4.4): the black-hole pair rate is the geometric middle of the gravitational-wave stock and the Lamberts simulation, which differ by 80 times. If you want the lower or the higher end, it is one number. Default: middle (0.2%).

## Evidence notes

[R] items to confirm when paper access is allowed: the Kelvin-Helmholtz time; the number of known high-mass X-ray binaries (about 130); the luminosity and mass of Cygnus X-1; the double pulsar's period and merger time; the abundance of barium stars among giants (about 1%); the planetary-nebula central star close-binary fraction (20%); that Hills ejection leaves a hypervelocity star single; the density of Milky Way-equivalent galaxies used to turn gravitational-wave rates into per-galaxy rates (0.0116 per Mpc^3, 1.16e7 per Gpc^3); the Lamberts et al. 2018 count of 1.2e6 black-hole pairs (read in a search-result summary of arXiv 1908.08775, not at source); the Gaia 0.1% prediction (a summary of arXiv 2502.03527); every Gaia and pulsar orbital figure.

Weakest numbers: the black-hole pair, black-hole with neutron star and black-hole with star fractions (a range of 80, 90 and 20 times); the direct-collapse mass that drives the wide-channel Monte Carlo; the 30% tertiary share; the compact-companion share of massive stars (0.3%). All are single entries in the dict, and a changed value changes the next scatter only.

[C] values come from scripts kept in the session scratchpad and not shipped: `counts.py`, `rows.py` (row counts and class mix from the repository's own `star_population` and `phenomenon_scatter` functions), `mc.py` and `mc2.py` (the survival Monte Carlo: Sana period law on 0.15 to 5.5 in log days, 90% of massive primaries with a companion, mass ratio uniform 0.05 to 1, thermal eccentricity above 7 days, instantaneous collapse on a random orbital phase with an isotropic Maxwellian kick, a system counted as interacting and excluded when its pericentre is under 9.3 AU, twice the 1,000 solar-radius giant phase), and `stock.py` (Peters axes and the stock of compact pairs from a merger rate: with a delay-time distribution proportional to 1/t between 50 Myr and 13.8 Gyr the unmerged tight pairs number the merger rate times 13.8 Gyr / ln(13.8 Gyr / 50 Myr) = 2.5 Gyr, and the share of a log-flat 0.003 to 100 AU draw that merges in a Hubble time is 0.19, 0.26 and 0.33). The Monte Carlo is a size check, not a population synthesis; it has no mass transfer, common envelope or tides.

## Sources

Seen in search-result summaries (no paper was opened): [Double neutron star merger rate, Pol et al. 2019](https://arxiv.org/pdf/1811.04086); [revisited 2020](https://arxiv.org/pdf/2107.13307); [formation of double neutron stars in the Milky Way](https://arxiv.org/pdf/2402.04658); [low neutron star merger rates](https://arxiv.org/pdf/2604.05059); [GWTC-4.0 population properties](https://arxiv.org/pdf/2508.18083); [Gaia NS1](https://arxiv.org/pdf/2402.06722); [Gaia black hole predictions](https://arxiv.org/pdf/2502.03527); [Gaia BH1 and BH2 origin](https://arxiv.org/pdf/2403.13579); [detached black hole-star binaries in clusters](https://arxiv.org/pdf/2306.13121); [detached black hole binaries with normal-star companions](https://arxiv.org/pdf/1909.11328); [stellar binaries that survive supernovae](https://arxiv.org/pdf/1810.08620); [supernova precursors with compact companions](https://arxiv.org/pdf/2401.02389); [binary black hole formation and metallicity](https://arxiv.org/html/2411.02484v1); [Moe and Di Stefano 2017](https://arxiv.org/pdf/1606.05347); [Offner et al. 2023](https://arxiv.org/pdf/2203.10066); [triple evolution (Toonen et al.)](https://ar5iv.labs.arxiv.org/html/1612.06172); [massive triples (Stegmann, Antonini and Moe)](https://arxiv.org/pdf/2112.10786); [black widows and redbacks](https://arxiv.org/pdf/2210.16295); [white dwarf binarity](https://arxiv.org/pdf/1209.1427); [Trapezium dynamical evolution](https://arxiv.org/pdf/1809.03537), [Orion Trapezium](https://arxiv.org/pdf/1701.03440); [hierarchical system architecture review](https://arxiv.org/pdf/2109.09118). Repository files read: `generation/phenomenon_scatter.py`, `bright_stars.py`, `star_population.py`, `run_sector.py`, `system.py`, `binary.py`, `wide_binary.py`, `tuning.py`; `docs/database-schema.md`, `docs/TODO.md` (GEN.129, GEN.130, GEN.223, GEN.250, GEN.251); the design notes linked above and `phenomenon-scatter-mass-cut.md`.
