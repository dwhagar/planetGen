# Size laws and satellite shares of small bodies (GEN.566, GEN.580)

Research note, no behaviour change. It lists each source with the numbers read in it and what is still missing, for the
belt member draw (GEN.566), the comet size cut of GEN.559 (GEN.560's interstellar size law), the hosting gate
`can_host(body)` (GEN.579) and the sub-satellite gate row (GEN.580). Raw papers stay out of the repo.

Source states: READ is a number found in the paper itself this session; 2ND-HAND is a number quoted by a paper that was read
(the original was not); RECALLED is from memory and not usable; STAND-IN is a placeholder named in the gap plan registry.
Rules (standing): verified wins; two verified in doubt, the newest wins; two in question, the widest range.

## 1. Size-frequency laws

Cumulative law used here: `N(>D) proportional to D^-q`. A magnitude law `dN/dH proportional to 10^(alpha H)` is the same law
with `q = 5 alpha` (Fraser et al. 2014 and Singer et al. 2019 both state the relation; a crater differential slope is
`q_diff = -(5 alpha + 1)`). Diameter from absolute magnitude: `D[km] = 1329 / sqrt(p) * 10^(-H/5)` (p the geometric
albedo).

| Population | Law | Range | Source | State |
|---|---|---|---|---|
| Main belt, small | q = 1.3 (differential index 2.3) | 0.4 to 5 km | Ivezic et al. 2001 (astro-ph/0105511), SDSS, 98 percent complete sample of about 10,000; corrected in GEN.1061, see 1.4 | READ |
| Main belt, medium | q = 3.0 (differential index 4.00 +- 0.05) | 5 to 40 km | Ivezic et al. 2001; corrected in GEN.1061, see 1.4 | READ |
| Main belt, exponents (GEN.1093) | checked verbatim in GEN.1061 (section 6.6 and Fig. 25 of the PDF): the 2.3 and 4 are indices of `dN/dD`, so the cumulative `q` is 1.3 and 3.0; the cumulative reading leaves 3 asteroids above 40 km against some hundreds in the belt | | Ivezic et al. 2001 | READ, section 6.6 |
| Main belt, count | about 530,000 objects with D > 1 km | belt | Ivezic et al. 2001 (four times fewer than earlier estimates) | READ |
| Main belt, small | q = 2.5, kink between 15 and 25 km | small side of the kink | Masiero et al. 2011 (arXiv:1109.4096); raw, not debiased (the paper says debiasing is needed) | READ |
| Jupiter-family comet nuclei | q = 1.92 +- 0.20 | r_N >= 1.25 km | Snodgrass et al. 2011 (arXiv:1101.4228); earlier values 1.45 to 2.7 depend on the size range chosen (Meech 2004: 1.45 over 1 to 10 km, 1.91 over 2 to 5 km: 2.7 for r >= 1.5 km) | READ |
| Cold classical KBOs | alpha1 = 1.5 (+0.4 -0.2), faint alpha2 about 0.2, break H_B = 6.9 (+0.1 -0.2) (r' band) | bright and faint sides | Fraser et al. 2014 (arXiv:1401.2157) | READ |
| Hot KBOs | alpha1 = 0.87 (+0.07 -0.2), faint alpha2 about 0.2, H_B = 7.7 (+1.0 -0.5) | | Fraser et al. 2014 | READ |
| Cold classical KBOs, newer | exponentially tapered or rolling power law fit equally well, single power law ruled out; M_CC(H_r < 12) = 0.0017 (+0.0010 -0.0004) Earth masses (p = 0.15, 1 g/cm3) | | DEEP survey paper V (arXiv:2309.09478) | READ |
| KBOs below about 2 km | crater differential slope q = -3 above 15 km craters, about -1.7 below; impactor magnitude slope alpha about 0.4 (+-0.04) above the elbow and about 0.15 (+-0.04) below, i.e. cumulative exponent about 2.0 above and about 0.75 below an elbow at D of 1 to 2 km | | Singer et al. 2019 (Science 363, 955; arXiv:1902.10795), New Horizons craters on Pluto and Charon | READ |
| Interstellar objects, 100 m scale | n = 0.1 per au3 (range 0.08 to 0.2) | about 100 m radius | Jewitt and Seligman 2023 (arXiv:2209.08182) quoting Jewitt 2017 (0.1) and Do et al. 2018 (0.2) | READ |
| Interstellar objects, 100 m scale | n = 0.2 per au3, mass density about 4 Earth masses per pc3 | | Do, Tucker and Tonry 2018 (arXiv:1801.02821) | READ |
| Interstellar dust to metres | flux compatible with r^-3 radius distribution | 0.1 micrometre up to metres | Jewitt and Seligman 2023 fig. 19 | READ |
| Interstellar objects, bodies | slope q constrained by detections only: q > 4 for asteroid-like, q > 4.5 for cometary at 0.1 per au3; TESS about 1 object for q < 2.5; LSST 0.2 to 100 per year over q 1 to 3 | | Albornoz-Montenegro and Fuentes 2026 (arXiv:2609.25488, 21 Sep 2026) | READ |
| Free-floating planets, mass | differential index p = 0.92 +- 0.06, size index 3.8 +- 0.2 | | Jewitt and Seligman 2023 | 2ND-HAND |

### 1.1 Computed from the table (checked, not read)

Fraser's breaks as diameters, `D = 1329 / sqrt(p) * 10^(-H/5)`:

| Population | H_B | p = 0.04 | p = 0.10 | p = 0.15 |
|---|---|---|---|---|
| cold | 6.9 | 277 km | 175 km | 143 km |
| hot | 7.7 | 192 km | 121 km | 99 km |

The cumulative exponents are `q = 5 alpha`: cold 7.5 bright side, hot 4.35, faint sides about 1.0; Singer's crater result
gives 2.0 above and 0.75 below the elbow. The interstellar density of 0.1 au^-3 is `0.1 / 1.1393e-16 pc^3 = 8.8e14 per
pc3` (range 0.08 to 0.2 au^-3 is 7.0e14 to 1.8e15 per pc3). Number above a 10 m size in the illustration of GEN.566 (a belt
of 1e6 bodies above 1 km): 1e10 at q = 2, 2.5e11 at q = 2.7, 1e12 at q = 3, so "1e10 to 1e12" in the item is right for
q 2 to 3 and about 2.5e11 at its stated 2.7.

### 1.2 Correction found while reading

`docs/design/unified-star-like-objects.md` line 140 states "the real density is about 1e12 per pc^3 (range 1e11 to 1e14)"
for interstellar objects, marked [S]. The measured 100 m scale density, 0.1 au^-3 (Jewitt and Seligman 2023) to 0.2 au^-3
(Do et al. 2018), is **7e14 to 1.8e15 per pc3**, three orders of magnitude higher. The conclusion of that paragraph (the
objects cannot be rows and the count below a lowest size is computed) holds with the higher number. TODO files this.

### 1.4 Correction: Ivezic's indices are differential (GEN.1061)

The first version of this table gave Ivezic et al.'s 2.3 and 4 as cumulative exponents `q`. The paper (section 6.6 and its Fig. 25) states them as the index `alpha` of `dN/dD proportional to D^-alpha`, which is a
cumulative `q = alpha - 1`: **1.3 below 5 km and 3.0 from 5 to 40 km**. (Taken as cumulative, 530,000 bodies above 1 km would leave only 3 above 40 km; as the paper has them it is about 130.) Found by
the Aux 7 lane in the tide work. The belt ledger (`physics/belt_ledger.py`) uses the corrected exponents; `q = 3` holds the same mass in every doubling of size, so the law has to steepen above 40 km, which the
measurement does not reach (stand-in BR20).

### 1.3 What is missing

- **Comets, long period and Oort cloud.** No size law was read for long-period comets (Francis 2005 and Fernandez and Sosa
  2012 are recalled and not read; the network did not return them from arXiv). GEN.559 needs this and waits on it.
- **Main belt below 400 m** and the near-Earth population as a whole (Bottke et al. 2005, Harris and D'Abramo 2015 are
  recalled, not read; neither is on arXiv).
- **Kuiper belt between the breaks.** Fraser 2014 gives alpha on each side; the diameter of the transition depends on the
  albedo, which is not one number (hot objects have a bimodal albedo, the paper says), and DEEP's tapered law is not given
  here as parameters. Read DEEP V's fit table before the draw uses it.
- **Dwarf and giant end.** The largest member of each population (Pluto-size) follows from the same laws: the hot law
  predicts a largest body about Pluto's mass (Bernstein et al. 2004, astro-ph/0308467, abstract READ).
- **Interstellar object size law** has no measurement, only detection constraints; the table above says what the surveys
  rule out. A draw cannot take one number as a result; it should take q from a range, widest first, 2.5 to 4.5, until a
  detection pins it. (This is a recommendation under the widest-range rule, not a read value.)

## 2. Binary and satellite fraction by size (GEN.580)

| Quantity | Value | Source | State |
|---|---|---|---|
| Binary fraction of NEAs above 300 m | 15 +- 4 percent | Pravec et al. 2006 via Margot et al. 2015 (Asteroids IV, arXiv:1504.00034) | 2ND-HAND (the review is READ) |
| NEAs above 200 m that are multiple systems | about one in six, from radar (Margot 2002, Taylor 2012) | Margot et al. 2015 | 2ND-HAND |
| Smallest radar primaries | 120 to 180 m (2003 SS84, 2004 FG11, 1994 CJ1) | Margot et al. 2015 | 2ND-HAND |
| Smallest photometric binary primary | 0.26 +- 0.03 km (2000 UG11); smaller ones exist but are hard to see | Margot et al. 2015 | 2ND-HAND |
| Largest photometric binary primary | about 13 km, (939) Isberga 13.4 +- 1.3 km | Margot et al. 2015 | 2ND-HAND |
| Main-belt binary fraction, small | not derived; discovery rate like the NEAs', so suspected similar; modelled steady-state 16 +- 6 percent | Margot 2015; Jacobson et al. 2016 (arXiv:1605.06157) | READ |
| Large main-belt asteroids with satellites | 16 known of over 300 surveyed; fraction "substantially smaller than the 16 percent" of NEAs; no satellite at any of 7 asteroids visited by spacecraft (Vesta, Lutetia, Steins, Toutatis, Annefrank, Itokawa, APL) | Margot et al. 2015 | 2ND-HAND |
| Jupiter Trojans with satellites | 2 confirmed (Patroclus, Hektor) | Margot et al. 2015 | 2ND-HAND |
| Satellite mass ratios, large main belt | 1e-6 to 1e-2 (Antiope excepted) | Margot et al. 2015 | 2ND-HAND |
| Cold classical KBOs with detectable binaries | 22 +10 -5 percent (inclination under 5 degrees, separations over 60 mas) | Noll et al. 2008 chapter (astro-ph/0703134) quoting Stephens and Noll 2006 | READ (primary 2ND-HAND) |
| Other dynamical classes combined | 5.5 +4 -2 percent (same limits) | same | READ |
| Known TNO binaries | 14 by adaptive optics, 64 by HST (2014) | Margot 2015 | 2ND-HAND |
| Origin of small binaries | near-critical angular momentum, rubble-pile primaries, rotational fission by YORP spin-up; separation lower limit a/D_p about 1.5 | Margot 2015; Jacobson et al. 2016 | READ |
| Fission mass ratio | high mass ratio (over 0.2) components four to eight times as frequent as low | Jacobson et al. 2016 | READ |

### 2.1 Sub-satellites of moons (the tidal gate row)

Kollmeier and Raymond 2019 (MNRAS, arXiv:1810.03304, READ):

- A submoon's orbit is stable out to a fraction f of its moon's Hill sphere; f = 0.4895 for prograde low-eccentricity orbits
  (Domingos et al. 2006, 2ND-HAND), about 0.4 when the planet's perturbation is included (Payne et al. 2013, 2ND-HAND).
- Tidal survival over T = 4.6 Gyr (their equation 1, adapted, read from the typeset page):
  `R_moon >= [ 39 M_sub k2 T sqrt(G) / (2 (4 pi rho_moon)^(8/3) Q) * (3 M_p / (f a_moon)^3)^(13/6) ]^(1/3)`, with the
  moon's Love number k2 (0.25), tidal quality factor Q (100), density rho_moon (2500 kg/m3), orbital radius a_moon, the
  planet mass M_p, the submoon mass M_sub (10 km radius, 2000 kg/m3 by default) and f the Hill fraction above. The
  dimensions work out (metres). Computed here with these defaults (not read): Moon 125 km, Callisto 331 km, Titan 354 km,
  Iapetus 35 km (all pass), Ganymede 1126 km (passes; the paper's text lists only the other four, so its figure uses
  something this note does not capture), Europa 3096 km, Io 8470 km, Rhea 2187 km and Triton 1489 km (all fail).
- Their result: 10 km submoons survive only around large (1000 km scale) moons on wide orbits. Among the Solar System's
  moons only Titan, Iapetus, Callisto and the Moon pass; none has a submoon, so formation or other effects (mass
  concentrations) also remove them. A mass ratio of about 1e-5 is a rough guide to long-lived submoons.

Missing: the tidal-stripping timescale by the planet for ordinary moons is only what that paper's equation gives; no
measured submoon exists to calibrate against.

### 2.2 Recommended rows for the hosting gate (GEN.579)

Under the widest-range rule, with the stand-in cut of Boss's card S5 (1 km) kept as the default:

| Row | Satellite share | State |
|---|---|---|
| Near-Earth-like rubble pile, D 0.2 to 13 km | 0.15 (range 0.11 to 0.19) | STAND-IN supported by a 2ND-HAND measure |
| Small main belt, D 0.3 to 13 km | 0.16 (range 0.10 to 0.22) | STAND-IN (modelled, not measured) |
| Large main belt, D over 100 km | 0.05 (range 0.02 to 0.16) | STAND-IN (16 of over 300, detectability-limited) |
| Cold classical KBO | 0.22 (range 0.17 to 0.32) | 22 +10 -5 percent, READ |
| Other trans-Neptunian classes | 0.055 (range 0.035 to 0.095) | READ |
| Below 0.2 km | 0.0 (not measured) | STAND-IN; evidence says satellites exist to about 120 m |
| Moon with submoon | 0.0 unless the tidal gate passes, then a stand-in of 1e-5 mass ratio | STAND-IN (Kollmeier and Raymond) |

The 1 km cut of card S5 is above every measured primary down to about 120 m, so it loses real satellites of small near-Earth
rubble piles; this is for the owner of GEN.579 to weigh (a row series by size class, as above, replaces the cut without
changing the gate's three parts).

## 3. Sources

READ this session (arXiv): astro-ph/0105511 (Ivezic 2001), 1109.4096 (Masiero 2011), 1101.4228 (Snodgrass 2011),
1401.2157 (Fraser 2014), 2309.09478 (DEEP V), 1902.10795 (Singer 2019), astro-ph/0308467 (Bernstein 2004, abstract),
2209.08182 (Jewitt and Seligman 2023), 1801.02821 (Do et al. 2018), 2609.25488 (Albornoz-Montenegro and Fuentes 2026),
1504.00034 (Margot et al. 2015), 1605.06157 (Jacobson et al. 2016), astro-ph/0703134 (Noll et al. chapter),
1810.03304 (Kollmeier and Raymond 2019).

Not read (recalled, not usable): Bottke et al. 2005, Harris and D'Abramo 2015, Francis 2005, Fernandez and Sosa 2012,
Pravec et al. 2006 and Margot et al. 2002 (cited only through the Asteroids IV review), Stephens and Noll 2006 (through the
chapter), Domingos et al. 2006 and Payne et al. 2013 (through Kollmeier and Raymond).
Masiero). None of these needs Boss; the first four are not on arXiv and were not fetched.
