# Captured asteroid moons: what the papers and the satellite tables say (GEN.1075)

Sources 1 lane, 10-30. This is the read list of [captured-asteroid-moons.md](../systems/captured-asteroid-moons.md) section 7 item 1: the papers opened, the head counts per giant read from the satellite tables, and what each changes in the design. The rows for GEN.725 are plan section 6.100 (CM1 to CM20). Docs only; no module reads these numbers yet.

## Result in one paragraph

The design's picture stands: windows only, the giants' four irregular populations are comparable, retrograde dominates 4:1, the outer edge is about half the Hill radius, and the surviving mass is about 1 percent of what was captured. Six numbers in the design change. (1) The capture probabilities of Saturn, Uranus and Neptune are not Hill-volume scalings; they are printed per planet and span a factor of 30 between instability histories, and the Jupiter table can be turned into a per-encounter rate the architecture can use (CM1, CM2). (2) The "63, 56, 27, 13" quoted as head counts are all moons, regular ones included; the irregulars were 55, 35, 9 and 8 in 2007 and are 107, 268, 10 and 8 today. (3) Above R = 5 km the giants hold 8, 11, 9 and 8 irregulars, so a Jupiter-like planet should come out near 8 (design 10 to 20). (4) Below 5 km the cumulative size law is 1.3 at Jupiter and 2.0 at Saturn, not the 2.5 the design had. (5) Storage above 1e16 kg (R = 10.6 km) is 5 to 8 moons per giant (design 8 to 12), and 90 percent of the mass is in 1 to 3 bodies. (6) The Kozai gap is hard at Jupiter and Saturn (0 of 375) and soft at the ice giants (4 of 18).

## Head counts per giant (satellite tables, read 10-30)

Counts are rows of the JPL planetary satellite mean-elements table (ssd.jpl.nasa.gov/sats/elem, 461 rows, provisional designations included); a moon is irregular when its mean semimajor axis is beyond the largest gap in `a` of its planet (Jupiter 3.0e6 km, Saturn 5.0e6, Uranus and Neptune 2.0e6 km; gaps of x3.1 to x15.5, so no moon sits near the cut; Triton and Iapetus are regular by this cut, Nereid irregular). Hill radius from the planet's mean distance and mass. Script and inputs: `/mnt/project-files/data/gen1075-captured-moons/analyse.py` (output beside it).

| | Jupiter | Saturn | Uranus | Neptune |
|---|---|---|---|---|
| Hill radius (million km) | 53.2 | 65.5 | 70.2 | 116.2 |
| all moons today (JPL) | 115 | 292 | 30 | 16 |
| irregular today | **107** | **268** | **10** | **8** |
| prograde / retrograde today | 15 / 92 | 58 / 210 | 1 / 9 | 4 / 4 |
| retrograde per prograde today | 6.1 | 3.6 | 9 | 1.0 |
| irregular in 2007 (JH07 Table 2: prograde + retrograde) | 55 (6 + 49) | 35 (8 + 27) | 9 (1 + 8) | 8 (4 + 4) |
| median a / r_H today (JH07 2007) | 0.43 (0.44) | 0.30 (0.29) | 0.15 (0.17) | 0.20 (0.19) |
| outer edge a / r_H, prograde / retrograde | 0.35 / 0.46 | 0.36 / 0.42 | 0.21 / 0.29 | 0.20 / 0.44 |
| irregular with inclination 55 to 130 deg | 0 | 0 | 1 (Margaret 60.5) | 3 (Halimede 119.6, Psamathe 127.8, Neso 128.4) |

Totals all four: 393 irregular, 315 retrograde and 78 prograde (4.0 to 1; JH07 88 to 19 = 4.6). JH07 Table 1 lists 55, 14, 9 and 7 irregulars; its own Table 2 gives 55, 35, 9 and 8, and the Neptune row of Table 2 does not add up (6 + 4 + 4 is not 13), so Table 2's irregular columns are used and the paper's Table 1 and Neptune total are taken as slips.

Size counts need the radii; these are Wikipedia's mean radii (from JPL and Sheppard's discovery papers; albedo 0.04 for most small ones), so they carry a factor-of-1.3 uncertainty in R. Complete to a radius of about 1.4 km at Jupiter, 1.6 km at Saturn, 7 km at Uranus ([R4053]: 50 percent at mR 26.1) and 16 km at Neptune (JH07: 36 km at mR 24, scaled to the 25.8 survey depth).

| irregulars with radius at least | Jupiter | Saturn | Uranus | Neptune |
|---|---|---|---|---|
| 1 km | 76 | 250 | 10 | 8 |
| 5 km | **8** | **11** | **9** | **8** |
| 10.6 km (1e16 kg at 2000 kg/m3) | **8** | **5** | **6** | **8** |
| 40 km | 1 (2 with Himalia 85 km and Elara 43 km, as JPL and NDD14 use) | 1 | 1 | 1 |
| largest | Himalia 70 | Phoebe 107 | Sycorax 93 | Nereid 178 (a / r_H 0.047, e 0.75, i 5 deg; without it Halimede 31) |
| total mass at 2000 kg/m3 (kg) | 3.9e18 | 1.0e19 | 7.5e18 | 4.8e19 (Nereid 98 percent) |
| share of the sum of R cubed in the largest | 0.74 | 0.98 | 0.90 | 0.98 |
| bodies holding 90 percent of the mass | 3 | 1 | 1 | 1 |
| cumulative index N(>R) ~ R^-b, 1 to 5 km | 1.3 | 2.0 | not complete | not complete |
| cumulative index 5 to 40 km | 0.9 | 1.3 | 1.0 | 1.0 |

## Papers on the read list

| Paper | Status | Taken |
|---|---|---|
| [R3265] (NVM07), AJ 133, 1962 | **read** (Troja copy) | per-planet capture probabilities: Uranus 2.6e-7, Neptune 5.4e-7 (means of 14 runs of the original Nice model), Jupiter 8.5e-9 and Saturn 2.4e-8 (job 9, the best of the runs; the mean at Saturn is 5.7 times lower); 35 earth masses of disc, size law q1 4.25, q2 2.75 differential, break D 100 km; largest captured D 135 km (Jupiter), 186 km (Saturn) against Himalia 170 and Phoebe 213 km; captured ice-giant populations 2 orders larger than the known ones |
| [R3264] (NDD14), [R3264] | **read** | Table 1 (Jupiter only; the design said per planet): 59, 280, 49 encounters, f3AU 4.2, 4.0, 3.5e-4, 1458, 3617, 1403 captured stable, Pcap 1.5, 3.6, 1.3e-8 (reproduced: Pcap = f3AU x Ncap / 4e7); Saturn about 5e-8, Uranus and Neptune (1 to 3)e-8 in the new histories |
| [R2158] (JH07), astro-ph/0703059 | **read** | Tables 1 to 3, retrograde to prograde 88 to 19, medians of a / r_H, no inclinations 55 to 130, Kozai limits 39.2 to 140.8 deg, Kozai times 180 yr (prograde) and 65 yr (retrograde) at Jupiter, family speeds about the escape speed of the largest member, collision with a regular moon per orbit (r / 2R)^2 |
| [R3147] (VNL08), AJ 136, 1463 | **read** (Troja copy) | binary exchange: Table 4; stable limits 0.4 r_H prograde and 0.5 to 0.6 retrograde (quoting [R3263]) |
| [R0483] | **read** | populations lose about 99 percent of their mass in hundreds of Myr and then sit in steady state; initial mass 0.0008 to 0.03 lunar masses (D above 0.1 km); Tables 1 and 2: collision probabilities and speeds; families |
| [R4053], astro-ph/0410059 (Uranus) | **read** | q = 1.8 +/- 0.3 differential above 7 km at Uranus, Jupiter 1.9 +/- 0.3 ([R4054] as quoted); completeness 7 km at Uranus |
| [R0163] (not on the list; the newest Jupiter size law) | **read** | luminosity index 0.29 +/- 0.15, q about 2.5 differential down to 0.4 km, about 600 (factor 2) retrograde moons above 0.4 km |
| [R0037] | **read** | capture of Triton from a binary of about 1000 km bodies; conditions (binary, disrupted, one bound) |
| [R3636] (not on the list; binary exchange in a migrating system) | **read** (abstract and introduction) | about 1 in 100 binary planetesimals gives a captured secondary at a Neptune-mass outer planet |
| [R2403] (not on the list; collisional capture in the gas era) | **read** (abstract and conclusions) | bound fraction of collisions in the Hill sphere at most 1e-5, 60 to over 1000 collisions per Myr, the colliders mostly destroyed |
| [R1666] | **read** (abstract) | temporary captures last 286 +/- 18 days and 2.88 +/- 0.82 revolutions; at least one 1 m natural Earth satellite at any time; 0.1 percent of meteors were captured; the one known case, 2006 RH120 |
| [R3814], Nature Geoscience (HAL copy) | **read** (abstract and introduction) | Phobos and Deimos are not captured: tidal dissipation cannot turn a captured body's orbit into their near-circular equatorial one; a 1e18 kg debris disc after a giant impact, inner large moons fall back in about 5 Myr |
| [R4456] (Saturn) | fetched, not used beyond JH07's families | Saturn groups of 3 and 4 prograde bodies |
| [R3263] | **not reached** (publisher copy refused) | second-hand: VNL08 (limits), JH07 (Kozai, a / r_H trend), [R0483] (families, collisions) |
| [R1771] | **not reached** | second-hand: JH07 section 5 (stability depends on e and i); the 0.4895 r_H circular limit is already in the repo |
| [R4053, R4054], 2006 | **not reached** (UCLA and publisher copies refused) | second-hand: JH07 Tables 2 and 3, [R0483, R4053] (Uranus survey read in full) |

## What this changes in the design

1. **Section 2.3 and the stand-in table, capture probability.** Replace "scaled by Hill volume to other planets" with the printed values (CM2). Ice giants capture more than gas giants in the original Nice model (3e-7 to 5e-7 against 1e-8 to 2e-8), and less in the newer histories ((1 to 3)e-8); the driver is the number of encounters the planet has with another giant, not its size. CM1 gives the per-encounter form (Pcap about 1.3e-10 to 2.6e-10 per encounter per disc particle at Jupiter, from 6e-7 stable captures per particle injected into the encounter sphere times the 3.5 to 4.2e-4 of the disc in it), so item 9's record of the encounter count of the architecture is enough to set a planet's capture probability.
2. **Section 2.3, head counts.** Use 55, 35, 9, 8 (2007) and 107, 268, 10, 8 (10-30) for irregulars; 63, 56, 27, 13 count regular moons too.
3. **Section 4.5, targets.** Irregulars above 5 km: Jupiter 8 (design 10 to 20), Saturn 11 (5 to 15), Uranus 9 and Neptune 8 (5 to 12; Neptune is complete only above 16 km, 7 there). Above 40 km: 1 at each (design 0 to 2). Largest: Himalia, Phoebe, Sycorax and Nereid, so the ice giant's largest is 93 km, or 178 km if Nereid is counted; the design's 30 to 100 km misses Nereid, which has the orbit of a scattered regular moon (a / r_H 0.047, i 5 deg).
4. **Section 2.3, size law.** Above 7 km the differential index is 1.8 to 1.9 (cumulative 0.8 to 1.0), as the design had. Below 5 km the design's "q about 3.5, cumulative 2.5" is not supported: counting the table gives 1.3 at Jupiter and 2.0 at Saturn, and the deepest Jupiter survey gives a differential 2.5 (cumulative 1.5) to 0.4 km. Use 1.3 to 1.5 at Jupiter and 2.0 at Saturn below 5 km; the 3.5 was marked "citation needed" in the source I had.
5. **Section 6, storage.** The floor of 1e16 kg (R = 10.6 km at 2000 kg/m3; 9.7 to 11.7 km for 2600 to 1500) stores 8, 5, 6 and 8 (design 8 to 12). Rule D's 90 percent of the mass is reached with 3, 1, 1 and 1 bodies, so the floor decides and the cap of 100 is never reached by the Solar System giants (it matters for a system with a heavy primordial reservoir).
6. **Section 2.1, Kozai band.** The edges 55 and 130 are the observed gap (JH07); the quadrupole limits are 39.2 and 140.8 degrees. Four of the 18 ice-giant irregulars lie inside the gap (Margaret at 60.5, and the three Neptune moons at 119.6 to 128.4 degrees, all at the far end of the retrograde zone); none of the 375 at Jupiter and Saturn do. `KOZAI_FORBIDDEN_DEG` is a hard cut at the gas giants and could be a soft one (a weight of about 0.2) at the ice giants.
7. **Sections 2.2 and 4.1, windows.** Gas-era collisional capture is negligible (CM12). Binary exchange gives under 0.001 small irregulars at Jupiter and Saturn, up to 0.07 at Uranus and up to 1.7 at Neptune if every planetesimal were a binary, and a chance of 0.2 to 3.5 percent of capturing a Triton-class body at Neptune, early (CM11). Temporary capture is a flag only (CM13).
8. **Collisions (item 5).** The present collisional lifetime of a 1 to 5 km irregular against all the others is 3e9 to 1.3e10 yr and of a 70 km one 9e7 to 1.3e9 yr (CM9), so the grinding is largely over today, and the 99 percent loss happened in the first hundreds of Myr [R0483]. The mass left now (4e18 to 1e19 kg at Jupiter, Saturn and Uranus) is 0.2 to 17 percent of the 0.0008 to 0.03 lunar masses (5.9e19 to 2.2e21 kg) the paper starts from, so it brackets the paper's 99 percent loss without fixing it. For the generator: draw the captured mass from the capture probability, then keep a fraction of 1 percent (0.5 to 5), with the largest body kept whole.
9. **Terrestrial planets.** [R3814] supports zero captures: the small circular equatorial moons of a rocky planet are impact debris, not captured bodies.
10. **Mass check of CM1.** (3 to 4)e7 disc particles above 40 km times 1.3e-8 to 3.6e-8 gives 0.4 to 1.4 captured bodies above 40 km at Jupiter (NDD14 says about one), against 1 or 2 seen.

## Not read, still open

- Gas-drag capture strength: see the GEN.1122 section below; the share stays a guess (CM14, burn-down line).
- Per-planet capture in the newer instability histories for Saturn, Uranus and Neptune beyond the one sentence of NDD14 section 3.4.
- Family speeds and counts at Saturn beyond JH07 and [R0483] (the Gallic, Inuit and Norse groups).
- The origin of Nereid (scattered regular moon after Triton's capture, or an irregular) decides whether Neptune's largest irregular is 178 km or 31 km.

## GEN.1122: gas-drag capture [R0956, R3560]

Sources 1 lane, 10-30. **Result: neither paper could be opened and no source reached gives a drag share against the encounter channel, so CM14 stays a guess (T6, L3); what the second-hand reading does give is a gate and a size window for the drag channel, now rows CM18 to CM20.**

What was reached: the abstract of [R3560] (Icarus 37, 587; NTRS, the full text is not on the server); [R0956] not on arXiv and its publisher page is refused, so read only through its two summaries: [R2158] section 6.2 (read in full) and [R4055], IAU 229, [R4052] section 3.1 (read in full). Also searched for a numeric capture efficiency in [R3513] (fetched, not read beyond its title: it is the three-body channel).

| Quantity | Value | Source |
|---|---|---|
| who can do it | gas giants only: Jupiter and Saturn kept about 300 and 100 earth masses of H2 and He, Uranus and Neptune about 1; the bloated young envelope was hundreds of times the planet's size | JH07 6.3 |
| capture condition | the body must meet about its own mass of gas on its way through the envelope | [R4055] 3.1 |
| size window | a few km to a few hundred km; larger bodies are hardly slowed, smaller ones spiral into the planet | [R3560] as quoted by [R4055] and JH07 6.2 |
| timing | the envelope must collapse within a few thousand years of the capture, else the body spirals in; the present irregulars would be the last few captured | [R4055] 3.1 |
| orbits it predicts | retrograde bodies lose more (smaller e, i, a), smaller bodies evolve faster and are removed first, larger ones are captured nearer the planet; none of these size trends is seen | [R4055] 3.1; JH07 6.2 |
| fit to the data | [R0956] explain the prograde Himalia group (disc, not envelope); JH07 call such a model dependent on weakly constrained parameters of the gas; the equal populations at the four giants argue against drag as the main agent | JH07 6.2, 6.3 |
| pull-down (sudden planet growth) | at least a 10 percent mass jump (Vieira [R4546]) or over 40 percent in a few thousand years; gas giants only; also works for prograde | JH07 6.1; [R4055] 3.2 |
| quasi-satellites | 5 to 20 percent of scattered planetesimals become quasi-satellites of the planet | JH07 6.3 |

For the generator: the drag channel is switched on only for a planet with a massive envelope in the gas era, for bodies of 1 to 500 km, and its yield stays the guess of CM14 (0.1 to 0.3 of the encounter channel) with the knob 0 to 0.3, since the observed equality of the four giants (JH07) and the lack of size-orbit trends [R4055] put the upper end at risk. A physical form is possible when the envelope is generated: capture if the column density of gas along the path times the cross-section is at least the body's mass, which gives R_max = (3/4) (rho_gas L / rho_body); it needs the gas density of the envelope (not in any paper reached).

Still open: the full text of [R0956, R3560] (a person with library access could read the capture-efficiency figure of the first and the disc-capture rate of the second).
