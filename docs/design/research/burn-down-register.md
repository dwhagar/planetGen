# The burn-down register: what each hand-tuned stand-in needs before it can go (GEN.946)

Status: hand-kept, checked by `src/tests/test_burn_down_register.py` and shown as the last column of `stand-in-values.md`. Every registry row at tweak level L3 or L4 (hand-tuned against a sweep, or fitted with no source) has a line here with one of three dispositions:

- `no data`: nothing we hold or can reach (primary or secondary) measures the quantity; the row keeps its tier and tweak level, and the reason says why in plain words. A reason is not a promise: the first measurement replaces the row.
- `replace`: a named open item replaces the number, because the data or the model it needs is that item's work. The item reads the data (Boss, 8:16 PM EDT 2026-10-29: secondary sources accepted) and lowers the row to L0 or L1.
- `replaced`: the number now comes from data; the line stays one release so the change is findable, and the row is no longer L3 or L4.

The test fails when a row at L3 or L4 has no line (so a new hand-tuned constant needs its disposition in the same change as its registry row, the same completeness rule as TEST.455 for unregistered constants), when a line names a row that is not at L3 or L4 and is not `replaced`, and when a `replace` item is not open in `docs/TODO.md`. The audit of every unregistered constant is GEN.945 (Bug 2); its high-ranked findings join this register as GEN.946 splits them into per-area items.

| row | disposition | item | reason |
|---|---|---|---|
| 5.2:BC6 | no data | - | exomoons are not observed; the regular satellites of the four giants are the only moons with a composition and the rule against fixed Solar System content keeps them a check; tier T6 stays |
| 5.3:AT.a1, a2 | no data | - | GEN.482 read the first-hand ranges at the Moon, Mars and Earth (T1 at those three points); only a larger white-dwarf ensemble improves the zero point |
| 5.3:AT.b5 | no data | - | Boss decided a log-uniform placeholder for sub-percent water on 2026-10-17; no ensemble of exoplanet water contents exists; any measurement replaces it |
| 6.7:AM64 | replace | GEN.804 | only figure read-offs of [R2184] exist; the upper-atmosphere model that reads the stored XUV replaces the fitted exobase constants |
| 6.7:AM65 | replace | GEN.804 | only figure read-offs of [R2184] exist; the upper-atmosphere model that reads the stored XUV replaces the fitted exobase constants |
| 6.7:AM66 | replace | GEN.804 | only figure read-offs of [R2184] exist; the upper-atmosphere model that reads the stored XUV replaces the fitted exobase constants |
| 6.7:AM67 | replace | GEN.804 | only figure read-offs of [R2184] exist; the upper-atmosphere model that reads the stored XUV replaces the fitted exobase constants |
| 6.7:AM71 | replace | GEN.804 | only figure read-offs of [R2184] exist; the upper-atmosphere model that reads the stored XUV replaces the fitted exobase constants |
| 6.7:AM81 | no data | - | the radius valley follow-up (GEN.722) is retired; the widths are fitted to the Kepler occurrence tables ([R1478, R4504] as second-hand widths) and any measured mass-loss width replaces them |
| 6.10:RV2 | no data | - | the radius valley follow-up (GEN.722) is retired; the widths are fitted to the Kepler occurrence tables ([R1478, R4504] as second-hand widths) and any measured mass-loss width replaces them |
| 6.10:RV3 | no data | - | the radius valley follow-up (GEN.722) is retired; the widths are fitted to the Kepler occurrence tables ([R1478, R4504] as second-hand widths) and any measured mass-loss width replaces them |
| 6.9:SP1 | no data | - | a display choice (how fast a rendered body turns on screen), not a physical quantity; tier T9 stays |
| 6.15:RB6 | no data | - | AP-9, AE-9 and the Divine-Garrett model are not held and the Galileo and Juno dose records are figures; the shapes stay recalled stand-ins |
| 6.15:RB7 | no data | - | AP-9, AE-9 and the Divine-Garrett model are not held and the Galileo and Juno dose records are figures; the shapes stay recalled stand-ins |
| 6.18:LF5 | no data | - | the pace timelines are scaled by hand from one history; the hard-step life timeline (GEN.750, closed) supersedes them for life and the table goes when the civilization check stops reading it |
| 6.18:LF6 | no data | - | design thresholds of the life scores (GEN.649): no measured quantity exists; tier T9 stays |
| 6.18:LF7 | no data | - | design thresholds of the life scores (GEN.649): no measured quantity exists; tier T9 stays |
| 6.18:LF8 | no data | - | design thresholds of the life scores (GEN.649): no measured quantity exists; tier T9 stays |
| 6.18:LF9 | no data | - | GEN.894 (merged, train 2) put the strange-solvent classes on real quantities, but habitability.SOLVENT_FACTOR is still read by the solvent gate; no measured limit of a non-water biochemistry exists |
| 6.18:LF24 | no data | - | GEN.894 merged (train 2) and the classes read real quantities; the strange-solvent classes replace the solvent factors and the methane-sea gate; no measured limit of a non-water biochemistry exists |
| 6.18:LF11 | no data | - | the habitability index's own scale choices (Boss, GEN.215, GEN.654); there is no data to read; tier stays |
| 6.18:LF13 | no data | - | the habitability index's own scale choices (Boss, GEN.215, GEN.654); there is no data to read; tier stays |
| 6.18:LF14 | no data | - | the habitability index's own scale choices (Boss, GEN.215, GEN.654); there is no data to read; tier stays |
| 6.18:LF26 | no data | - | the habitability index's own scale choices (Boss, GEN.215, GEN.654); there is no data to read; tier stays |
| 6.18:LF28 | no data | - | the habitability index's own scale choices (Boss, GEN.215, GEN.654); there is no data to read; tier stays |
| 6.18:LF15 | no data | - | phosphate classes of an ocean: Earth's ocean is the only measured one; the other classes are chosen |
| 6.18:LF16 | no data | - | phosphate classes of an ocean: Earth's ocean is the only measured one; the other classes are chosen |
| 6.18:LF17 | replace | GEN.878 | a curve per metabolism and the energy source per class (pigment band times the star's spectrum) replace the one log scale and the shares |
| 6.18:LF18 | replace | GEN.878 | a curve per metabolism and the energy source per class (pigment band times the star's spectrum) replace the one log scale and the shares |
| 6.18:LF23 | no data | - | the animal oxygen floor by grade is a literature number ([R3084], recalled) for one grade; GEN.877 is closed and kept the band; the widths are chosen |
| 6.18:LF25 | no data | - | one score per rung of Table 2 of the Exoplanetary Habitability doc; design numbers with no data |
| 6.18:LF49 | no data | - | the hard-step means follow [R4642] in form; the only calibration is one planet's history, and the rule against fixed Solar System content keeps it a check |
| 6.18:LF50 | no data | - | the hard-step means follow [R4642] in form; the only calibration is one planet's history, and the rule against fixed Solar System content keeps it a check |
| 6.18:LF51 | no data | - | the hard-step means follow [R4642] in form; the only calibration is one planet's history, and the rule against fixed Solar System content keeps it a check |
| 6.18:LF57 | no data | - | fitted to the first census (TEST.430) by design: the prevalence targets are the only calibration, so the rows stay L3 until a measured fraction of worlds with each rung exists |
| 6.18:LF58 | no data | - | fitted to the first census (TEST.430) by design: the prevalence targets are the only calibration, so the rows stay L3 until a measured fraction of worlds with each rung exists |
| 6.19:LF46 | no data | - | fitted to the first census (TEST.430) by design: the prevalence targets are the only calibration, so the rows stay L3 until a measured fraction of worlds with each rung exists |
| 6.19:LF31 | no data | - | the pigment light windows are published absorption bands; the gates and the temperature windows are chosen limits of the one known biosphere |
| 6.19:LF32 | no data | - | the pigment light windows are published absorption bands; the gates and the temperature windows are chosen limits of the one known biosphere |
| 6.19:LF33 | no data | - | the pigment light windows are published absorption bands; the gates and the temperature windows are chosen limits of the one known biosphere |
| 6.19:LF34 | replace | GEN.878 | energy by chemistry and the stored supplies (GEN.887, closed) replace the placeholder supplies of the classes; no measured flux of these metabolisms on other worlds exists |
| 6.19:LF35 | replace | GEN.878 | energy by chemistry and the stored supplies (GEN.887, closed) replace the placeholder supplies of the classes; no measured flux of these metabolisms on other worlds exists |
| 6.19:LF36 | replace | GEN.878 | energy by chemistry and the stored supplies (GEN.887, closed) replace the placeholder supplies of the classes; no measured flux of these metabolisms on other worlds exists |
| 6.19:LF37 | replace | GEN.878 | energy by chemistry and the stored supplies (GEN.887, closed) replace the placeholder supplies of the classes; no measured flux of these metabolisms on other worlds exists |
| 6.19:LF38 | replace | GEN.878 | energy by chemistry and the stored supplies (GEN.887, closed) replace the placeholder supplies of the classes; no measured flux of these metabolisms on other worlds exists |
| 6.19:LF39 | replace | GEN.878 | energy by chemistry and the stored supplies (GEN.887, closed) replace the placeholder supplies of the classes; no measured flux of these metabolisms on other worlds exists |
| 6.19:LF40 | no data | - | GEN.894 merged (train 2) and the classes read real quantities; the strange-solvent classes 5 to 9 carry placeholder supplies and windows; there is no organism of these classes to measure |
| 6.19:LF41 | no data | - | GEN.894 merged (train 2) and the classes read real quantities; the strange-solvent classes 5 to 9 carry placeholder supplies and windows; there is no organism of these classes to measure |
| 6.19:LF42 | no data | - | GEN.894 merged (train 2) and the classes read real quantities; the strange-solvent classes 5 to 9 carry placeholder supplies and windows; there is no organism of these classes to measure |
| 6.19:LF43 | no data | - | GEN.894 merged (train 2) and the classes read real quantities; the strange-solvent classes 5 to 9 carry placeholder supplies and windows; there is no organism of these classes to measure |
| 6.19:LF44 | no data | - | GEN.894 merged (train 2) and the classes read real quantities; the strange-solvent classes 5 to 9 carry placeholder supplies and windows; there is no organism of these classes to measure |
| 6.19:LF45 | no data | - | GEN.894 merged (train 2) and the classes read real quantities; the strange-solvent classes 5 to 9 carry placeholder supplies and windows; there is no organism of these classes to measure |
| 6.19:LF47 | no data | - | chosen or Earth-calibrated numbers of the life supplies and ledgers (each carries its own tier in its file); no other biosphere exists to calibrate them |
| 6.19:LF59 | no data | - | chosen or Earth-calibrated numbers of the life supplies and ledgers (each carries its own tier in its file); no other biosphere exists to calibrate them |
| 6.19:LF60 | no data | - | chosen or Earth-calibrated numbers of the life supplies and ledgers (each carries its own tier in its file); no other biosphere exists to calibrate them |
| 6.19:LF61 | no data | - | the gates and the liquid window of the acid cloud are chosen limits and the acid strength line is fitted through two Venus cloud levels; no second acid world exists to calibrate against, and the [R4877] curve is not in the reachable papers |
| 6.19:LF62 | no data | - | the ammonia share of an icy body's ices is a choice between the comet value and the cosmic cap; no survey measures the ammonia of the ocean of any world but a few plume spectra, and the Mg/Si line is the textbook one |
| 6.19:LF53 | no data | - | chosen or Earth-calibrated numbers of the life supplies and ledgers (each carries its own tier in its file); no other biosphere exists to calibrate them |
| 6.20:PM2 | no data | - | the observed occurrence rates [R1156, R1479, R3499] are the bands the planet count and spacing fit to (GEN.728, GEN.729, GEN.828, GEN.829, built); the steepening itself is tuned to them and no disc survey measures it |
| 6.21:BO1 | no data | - | the oxygen balance has one measured world; Earth's fluxes are recalled [R1859, R1980] and calibrate the scale; tier stays |
| 6.21:BO2 | no data | - | the oxygen balance has one measured world; Earth's fluxes are recalled [R1859, R1980] and calibrate the scale; tier stays |
| 6.21:BO3 | no data | - | the oxygen balance has one measured world; Earth's fluxes are recalled [R1859, R1980] and calibrate the scale; tier stays |
| 6.21:BO4 | no data | - | the oxygen balance has one measured world; Earth's fluxes are recalled [R1859, R1980] and calibrate the scale; tier stays |
| 6.21:BO6 | no data | - | the oxygen balance has one measured world; Earth's fluxes are recalled [R1859, R1980] and calibrate the scale; tier stays |
| 6.21:BO8 | no data | - | the oxygen balance has one measured world; Earth's fluxes are recalled [R1859, R1980] and calibrate the scale; tier stays |
| 6.22:CC1 | no data | - | GEN.962 built the moon count from giant-impact frequency (rows RM1 to RM3), but no exomoon mass distribution is observed; the mass floor and the mass function stay recalled, with the satellite mass ratios a Solar System check, not a source |
| 6.22:CC2 | no data | - | GEN.962 built the moon count from giant-impact frequency (rows RM1 to RM3), but no exomoon mass distribution is observed; the mass floor and the mass function stay recalled, with the satellite mass ratios a Solar System check, not a source |
| 6.22:CC4 | replace | GEN.684 | albedo and greenhouse from the composition (the ground and cloud albedo) delete the regime ranges |
| 6.22:CC6 | no data | - | the lightest planet is a naming boundary and the habitable-request steering is a user-interface choice; neither is a measured quantity |
| 6.22:CC10 | no data | - | the lightest planet is a naming boundary and the habitable-request steering is a user-interface choice; neither is a measured quantity |
| 6.23:TR3 | no data | - | the form follows from physics (a suspension has no static rigidity); the floor is the lowest value that keeps the Io state's k2 under the fluid limit, and no measurement of the rigidity of a partial melt at planetary scale exists |
| 6.23:TR10 | no data | - | measured tidal Q of cold rock (Moon 38 to 41, Mars 86) lies far under the ceiling, so it constrains the ceiling only from below; the seismic Q the row cited is another frequency; no tidal Q above 100 is measured; tier T6, tweak L3 stay |
| 6.26:GS4 | no data | - | the distribution of the last impactor by planet mass is in [R0036, R2432] as figures and fits we do not hold as a table; a table of the last impact's mass ratio by planet mass replaces the isolation mass and steepness (the planet-formation train can give the isolation mass of the disc, formation.isolation_mass_kg) |
| 6.26:GS5 | no data | - | the same N-body table replaces the log-uniform range of the last giant impact |
| 6.26:GS6 | no data | - | the same N-body table (the planetesimal-built embryos) replaces the log-uniform range |
| 6.36:PA1 | no data | - | no model of the outer disc solids is filed; the numbers are tuned to the giants at 2 to 8 AU [R1479], and a measured pebble-flux model replaces them |
| 6.36:PA2 | no data | - | the share of cores that find gas is tuned to the giant occurrence and its rise with [Fe/H] [R2172]; no survey measures the disc lifetime it stands for |
| 6.36:PA3 | no data | - | the starting count of the inner chain is tuned to the [R4898] multiplicities after the formation-era dynamics build (GEN.730); the starting count of a real disc is not observed |
| 6.36:PA4 | no data | - | the masses of the Kepler planets are measured for few; the scale is tuned to the 1 to 10 Earth-mass range of [R4897] and its trend with host mass and metallicity |
| 6.36:PA5 | no data | - | [R4660] measure the pair ratio, not the spread of masses a system starts with; tuned to the 65 percent outer-larger share |
| 6.36:PA6 | no data | - | the starting spacing is tuned to the Pu and Wu 2015 and [R4660] spacings after the instability build (GEN.730); the spacing before the instability is not observed |
| 6.36:PA7 | no data | - | the inner edge of the chain is measured as a fall of occurrence inside 10 days [R3196]; the median period and width are tuned to the Kepler rise with period |
| 6.36:PA9 | no data | - | the share 0.0092 now sits between two read statistics ([R4897]: 0.62 percent Kepler, 0.9 to 1.2 percent radial velocity), but the row also holds the metallicity exponent 2.2 and the period range, which no survey constrains separately; the row stays L3 until they are split |
| 6.36:PA10 | no data | - | the inward move of other giants is tuned to the 1 to 3 AU peak of [R0960, R1479]; no survey measures the migration |
| 6.36:PA11 | no data | - | [R4271] (read) give the planet counts per log mass ratio beyond the snow line down to their sensitivity edge but no floor below it; the lightest oligarch that counts stays tuned to those counts; tier T5, tweak L3 stay |
| 6.36:PA14 | no data | - | a requested planet class takes its mass log-uniform over 0.05 to 1000 Earth masses after 8 free draws; no survey measures the mass of a planet picked by class, so the range is a design choice kept wide |
| 6.40:RD1 | no data | - | no atmosphere of a runaway-greenhouse planet has been measured; TRAPPIST-1 b and c only forbid keeping all of it |
| 6.40:RD2 | no data | - | no abiotic oxygen air has been detected or ruled out on an M-dwarf habitable-zone planet; Luger and Barnes fix only the possible amount |
| 6.40:RD3 | no data | - | no measurement of rocky exoplanet pressures by depth below the shoreline exists; GEN.483 derives the caps at these edges from the Solar System |
| 6.40:RD5 | no data | - | the 100 K is the width over which the equilibrium lines move with pressure and metallicity; a population of hot-giant chemistry regimes has not been measured |
| 6.41:TK2 | no data | - | the line between hydrated and primitive bodies is fitted to the C to P change of the belt; no model of aqueous alteration against ice content reaches it, and no belt outside the Solar System is resolved into taxa |
| 6.41:TK5 | no data | - | the accretion spread was recalled and is now read: [R2677] give core formation under about 1 Myr in the first reservoir and 1.8 to 4 Myr in the second, consistent with the 3 Myr; the value remains the one Solar System spread tuned by the melting test, no other system dates its planetesimals; tier T3, tweak L3 stay |
| 6.43:BM2 | no data | - | the words dense, sparse and typical are labels for the spread of belt masses at one age; no survey resolves the masses of belts of asteroids, and the factor of 3 is a choice of where the words change |
| 6.44:MG2 | no data | - | the widths of the axes and the distance past which a body is called none are choices; the collection's spread inside each group is not in the table, which gives averages and ranges only |
| 6.47:SE4 | replaced | - | the tail now comes from the isotope record of [R4475] (0.0077 events a year at 1e10 cm-2, power law -4 beyond); no longer hand-tuned; the pairing of fluence with peak flux stays McCracken's recalled number (tier T3, tweak L1) |
| 6.14:SP22 | no data | - | the generator holds no secular frequencies or fossil shapes to give the Cassini-state share of the inclination; a Laplace-Lagrange run of the generated system with the body's gravity field replaces the log-uniform share, which is kept inside Mercury's 0.004 and the bound 1
| 6.51:PS14 | no data | - | the label is a margin between the model's own solar-metallicity and metal-poor stars; no survey count of subdwarfs by luminosity deficit is held, and the first one replaces the margin |
| 6.60:GM3 | no data | - | a design default of the multi-star build; no measured limit for a close pair's separation in summed radii is held, and an N-body audit (GEN.1006 section 9) is optional; tier T6 stays |
| 6.60:GM7 | replace | GEN.1006 | the 3+1 network of [R4575] replaces the margin inside its training range; outside it the margin stays, flagged lower confidence |
| 6.60:GM12 | replace | GEN.1016 | the criterion holds only for q_out up to 5; GEN.1028 enforces the range (built), and GEN.1016 replaces it with the quadruple fit and the 3+1 criterion beyond it |
| 6.60:GM13 | replace | GEN.1006 | the 3+1 axis margin of 1.3 is hand-tuned on 331 N-body systems; a larger own N-body sample (GEN.1006 section 9) or a retrained network replaces it |
| 6.60:GM14 | replaced | - | GEN.1155: the 2+2 axis margin is now 1.5, tuned (L2) to the paper's 0.94 precision on 302 + 185 own REBOUND 2+2 systems; it was the 3+1 margin of 1.3 with no 2+2 sample |
| 6.148:RGD1 | no data | - | the clean-ice reflectance and the half-darkening share of an ice and dark mixture are recalled ([R0814] was not read); the share is tuned to Saturn's ring reflectances, and a measured mixture curve replaces both |
| 6.71:RR1 | no data | - | the rate of a resonance ring on a small body has no survey behind it: the three known rings are of bodies chosen for occultation campaigns, and the count of bodies searched is not published as a rate |
| 6.71:RR4 | no data | - | the mass and width range of a thin resonance ring is bracketed by the optical depth it must give; the measured optical depth and width of Quaoar's ring were not read |
| 6.69:LQ3 | no data | - | the brine and ammonia-water eutectic densities, viscosities and boiling points are recalled values; laboratory brine tables (Pitzer-model compilations) can replace them when read |
| 6.69:LQ5 | no data | - | the molten iron sulfide and Fe-S eutectic melt properties are recalled; measured Fe-S melt properties (Terasaki et al.) replace them when read |
| 5.2:BC20 | no data | - | the share of giants with a ring is now an output of the fall and passage histories (`ring_origin`); no ring census of exoplanets is held, and Cassini's ring age and mass are argued, so the history's own stand-ins RO1 to RO6 carry it |
| 5.2:BC21 | no data | - | the edges of ring systems are measured for four giants only; the outer edge against the Roche limit and the inner edge in planet radii have no second sample |
| 5.2:BC22 | no data | - | the kind limits are the middles of the gaps between four Solar System rings; a second sample of measured ring masses replaces them |
| 5.2:BC23 | no data | - | the share a ring begins with is still drawn from this range (no measured ice share of a young ring exists); what the dust has done to it since is worked out (`rings.darken`, GEN.1055 R6) |
| 6.70:RO1 | no data | - | no simulation table of the ring left by a disrupted moon by composition is in hand; [R0414, R0637] disagree on the share |
| 6.70:RO2 | no data | - | the age of Saturn's rings is contested between 1e8 years and the age of the planet, and no other ring system is dated |
| 6.70:RO3 | no data | - | no survey measures how often a ring-making comet passes inside a giant's Roche limit |
| 6.70:RO5 | no data | - | no N-body table of the share of a torn comet that stays in orbit by speed and periapsis is in hand; [R1140] gives a range only |
| 6.70:RO13 | no data | - | no impact-ejecta experiment on ring ice at ring meteoroid speeds gives the mass yield; the ejecta mass is argued over two decades |
| 6.70:RO14 | no data | - | no survey of the inner moonlets of giants outside the four of the Solar System is in hand |
| 6.70:RO16 | no data | - | no impact-ejecta experiment on a small moon's surface gives the kept share of the dust |
| 6.74:SH3 | no data | - | the pressure a given shoreline ratio allows is not measured as a table: only Mars's 6 mbar and Pluto's 10 micro-bar sit at ratios of 1.8 and 0.78; the caps stay the design note's steps until a survey of thin atmospheres fits them |
| 6.78:GT3 | replace | GEN.1058 | the share of ice giants with a recorded merger and the impactor mass share are tuned to one Solar System count (Uranus); the impact history of the outer chain, once recorded, gives the merger and its impactor |
| 6.96:IO2 | no data | - | fitted with the Jupiter row to the retrograde to prograde split; a capture simulation's width replaces it when one is held |
| 6.96:IO4 | no data | - | no family velocity field is held beyond the recalled tens to hundreds of m/s; the first one replaces the range |
| 6.83:SX2 | no data | - | no simulation or survey gives the captured share of a swept span for a given migration; replaced by the resonant populations of a Neptune-migration run |
| 6.83:SX3 | no data | - | the ejected share of a scattering giant is not tabulated by Safronov number; replaced when the scattering simulations are read |
| 6.83:SX4 | no data | - | no survey statistic gives a hot to cold share for an arbitrary system; the observed classical populations are the only count |
| 6.83:SX5 | no data | - | the split of the bound swept mass among scattered, centaur, detached, trojan, implanted and moon populations is not measured; replaced by a population synthesis |
| 6.95:BR20 | no data | - | the size-frequency of the main belt above 40 km ([R0484] and the catalogue counts) was not reachable: the paper is not on arXiv and no catalogue was fetched; Ivezic et al. measured only to 40 km. Any census of bodies above 100 km replaces the break and the steep slope |
| 6.100:CM14 | no data | - | the share a gas disc adds by drag is in no paper reached: [R3560] is the abstract only (NTRS), [R0956] is not reachable, and the two reviews read [R2158, R4055] give the gate and size window (CM18) and call the model weakly constrained; the first measured capture simulation of the drag channel replaces it |
| 6.98:BW2 | no data | - | no measurement of surface winds across the thin-air range (Mars, Pluto, Triton) gives the fade's shape; the first drag-law fit from surface wind data replaces it |
| 6.98:BW4 | no data | - | tuned to the Earth's mid-latitude surface wind; a boundary-layer drag law on the surface roughness replaces it when the surface roughness is stored |
| 6.98:BW5 | no data | - | the Earth's night to day wind ratio is the only value held; a stable boundary-layer model replaces it |
| 6.87:RK7 | no data | - | the stored Loss record holds no closest approach (distance, formation distance, time, scatterer orbit and mass, Safronov number, Hill radius only); the curve through three measured points stays until a deterministic closest approach is derived from the scatterer geometry |
| 6.112:MO2 | no data | - | the stored Loss record holds no closest approach; the curve through two printed points stays until a deterministic closest approach is derived from the scatterer geometry |
| 6.86:TU1 | no data | - | the level ceilings are the numbers of the Rogues dropdown that Boss sets, a display choice and not a research rate; tier stays |
| 6.86:TU2 | no data | - | uniform placement is the neutral choice; the smoothed rogue field (GEN.1115) is the only thing that could change it |
| 6.86:TU3 | no data | - | top-up rogues have no donor and no ejection record by design (apparent origin), so there is no scatterer distance to read; the drawn range stays |
| 6.108:TF1 | replace | GEN.1103 | a physical source for the mass the Rogues level cannot hold (the unbound remnants of the channels of section 6.103) replaces the stand-in origin |
| 6.108:TF2 | replace | GEN.1105 | the size-frequency law of the collisional fields sets the largest body and the density of the pieces |
| 6.108:TF3 | no data | - | nothing measures the spread of a field that was never observed; the dispersion of a real source (IF12) is the only thing that could replace it |
| 6.108:TF4 | no data | - | the counts and masses are chosen to keep a sector's top-up to a dozen clouds; no measurement exists to replace them |
| 6.97:BG14 | no data | - | no survey gives the mass of the warm dust of an inner belt well enough to set the least mass a belt keeps; the value is set so the share of main-sequence systems with an inner belt lands near 9 percent, inside the 4 to 14 percent of GEN.1156, and a size law with each belt's own largest body (GEN.566) or a dust survey would replace it |
| 6.117:BHS3 | no data | - | no source gives the share of a belt (as against a planet system) that a star's sudden mass loss destabilises; the value is set so the white dwarf pollution the model gives is in the observed range, and any N-body run of a belt round a mass-losing star replaces it |
| 6.118:IP5 | no data | - | no paper gives the gas-drag share against the three-body channel ([R3560] only as an abstract, [R0956] not reached); a gas-drag capture run of the generated envelope replaces it |
| 6.118:IP8 | no data | - | no binary-exchange run gives the size of the captured body for a general disc; the design's few hundred km is the only guide |
| 6.118:IP10 | no data | - | the disruption size ratio is an energy argument from recalled strengths, not a measurement; the size-dependent strength of a cascade run of the generator's own population replaces it |
| 6.119:CP3 | no data | - | no source gives how many pieces a captured parent was made of; the Himalia group's four known members are the only guide, and a parent built from the disc's own bodies (belt design) replaces it |
| 6.121:TJ4 | no data | - | no source read gives the leading share of a swarm; the capture history of the generated giant replaces it |
| 6.121:TJ5 | no data | - | no source read gives the eccentricity and inclination limits; an N-body map of the trojan region by mass ratio replaces them |
| 6.41:TK6 | no data | - | the within-family shares are fixed fractions; the meteorite-group map of GEN.475 labels the group of a mix but gives no shares for the taxa of a family, and no table of asteroid taxa by mass in a belt exists to fit them to |
| 6.130:AB1 | replace | GEN.804 | the exponents are tuned to the earlier table's level; the rotation and X-ray evolution of [R2183] (the model, not its figure) gives the decay itself, and the upper-atmosphere model that reads the stored XUV is the item that needs it |
| 6.43:LB1 | replaced | - | GEN.1000: each belt's largest body now comes from the size law (`belt_ledger.belt_state`: drawn from the belt's own stream and settled with the collisional mass); the range stand-in is gone |
| 6.45:MF1 | no data | - | no measurement gives the melt share of a mantle outside the Earth (and Io's lavas); a mantle potential temperature and solidus model is the replacement and nothing we hold supplies it |
| 6.125:TB3 | no data | - | no source read gives the share of co-orbital solids that grow into a trojan; a disc formation map replaces it |
| 6.125:TB4 | no data | - | no source read gives the shortest disc life for a trojan to grow; a formation map replaces it |
| 6.150:CF2 | no data | - | the 0.5 g and 2 g lines of the habitable classes have no measured soft edge; the width is hand-picked until the GEN.1207 sweep shows how far the class scatter reaches |
| 6.150:CF3 | no data | - | the temperature shoulder is hand-picked; a measured scatter of one class at fixed insolation would replace it |
| 6.150:CF4 | no data | - | the pressure shoulder is hand-picked on a log scale; no source gives a soft edge for the survival pressure limits |
| 6.150:CF5 | no data | - | the ocean-fraction shoulder is hand-picked; the law of GEN.1147 fixes how bimodal the cover is and the GEN.1207 sweep then sets the width |
| 6.150:CF6 | no data | - | the airless icy line shoulder is hand-picked; no source gives a soft edge for dry against icy airless bodies |
| 6.150:CF7 | no data | - | the gas-share shoulders are hand-picked; no source gives a soft edge for the oxygen, carbon dioxide and reducing-gas limits |
| 6.155:LR1 | no data | - | the abundance word edges on the 0 to 9 kind score are a hand-picked reading of a standard ecology scale (DAFOR style); no source ties the words to our score |
| 6.155:LR2 | no data | - | the trend look-ahead, band and spread are hand-picked; no source gives a rate of change that counts as a trend for a living class |
| 6.155:LR3 | no data | - | the time-left limits for Dying out, Precarious and Marginal are hand-picked; no source gives them for a living class |
| 6.155:LR5 | no data | - | the milestone each kind spreads after is a hand-picked reading of the life timeline; the stored schedule of the world replaces it where there is one |
| 6.126:RN3 | no data | - | no measured strength of ice with rock grains above 56 percent solids by volume; a measured curve replaces the extrapolated leg from 0.56 to the rock strength at 1 |
| 6.126:RN5 | no data | - | no source gives a porosity law for body-scale strength; the rubble-pile branch (friction plus a little cohesion) replaces it |
| 6.129:IT9 | no data | - | no measurement of the ammonia share at which a water ice mix takes the dihydrate's strength |
| 6.59:CK1 | replace | GEN.1021 | the N2 and Ar trapping shares are tuned to the one comet measured (67P); a second comet's N2 or Ar per water replaces them |
| 6.59:CK2 | replace | GEN.1021 | the OCS and SO2 trapping shares are tuned to the one comet measured (67P); a second comet's sulfur species replace them |
| 6.59:CK5 | replace | GEN.1021 | O2, krypton and xenon per water are the 67P values with no host ratio; a second comet's values replace them |
| 6.136:CC7 | no data | - | no measured viscosity of the lowest mantle at the core-mantle boundary; the Earth's core temperature and inner core size are the only constraint it is tuned to |
| 6.136:CC9 | no data | - | no measurement of the shell thickness at which a core dynamo dies; Mars and the Earth are the only two points |
| 6.142:DY7 | no data | - | the core model has the pure-iron liquidus only; no sulfur-bearing melting curve is held, and the first one replaces the stand-in moment of a small iron-rich body |
| 6.160:LUM2 | no data | - | the roundness a lump is drawn with when none is stored is a picture choice; the stored ratio of every system made after GEN.1159 replaces it, and nothing measures the shape of a body that was never stored |
| 6.160:LUM4 | no data | - | the bump height is tuned by eye; a stored relief or a fitted power spectrum of small-body topography (no source in hand) would replace it |
| 6.161:FM1 | no data | - | the line is the belts' (BG14): no survey measures the dust of a field; it goes where BG14 goes |
| 6.161:FM2 | no data | - | no model or measurement of the size distribution of the rubble of a planet-planet collision is in hand; the rocky main-belt law stands in until the largest-fragment reading (LS12) |
| 6.166:IHA1 | no data | - | no measured viscosity of any rocky planet's mantle other than the Earth's is in hand; the prefactor is retuned to hold the reference viscosity, not fitted to a planet |
| 6.53:DU1 | no data | - | no survey measures the dust of the disks this generator builds; the limit at the Sun's temperature stands where the A, F, G and K dwarfs hold a disk as often as the surveys find (GEN.1278 made it fall with the cube of the star's temperature and summed the dust of a wide belt's rings), and goes when the outer disk's mass and width come from a planetesimal-formation model or each survey's own flux limit is used |
| 6.185:BPO1 | no data | - | no survey of another star's outer disk gives the orbits of its populations; the cold stirring stands until a self-stirring or capture model gives it |
| 6.185:BPO2 | no data | - | no survey of another star's outer disk gives the orbits of its populations; the resonant eccentricity range stands until a self-stirring or capture model gives it |
| 6.185:BPO3 | no data | - | no survey of another star's outer disk gives the orbits of its populations; the crossing perihelion range stands until a self-stirring or capture model gives it |
| 6.185:BPO4 | no data | - | no survey of another star's outer disk gives the orbits of its populations; the detached perihelion range stands until a self-stirring or capture model gives it |
| 6.183:DM1 | no data | - | no survey gives the climate scatter among bodies of one class; the 5 K shoulder of the regional temperature lines is a width of the class fit, chosen so a class edge is not a cliff |
| 6.183:DM2 | no data | - | no measurement of another world's wet-bulb limit is in hand; the 3 C shoulder is a width of the class fit around the human limit of 31 C |
| 6.183:DM3 | no data | - | no class reads a humidity edge yet and no survey measures one; the shoulder waits for the first class that does |
| 6.183:DM5 | no data | - | the edges of classes M, H, O, P and R are the plan's reading of the class definitions; no survey of other worlds fits them, so they stay as fitted to the finished bodies (reports/class-fit-shadow-2026-10-31.md) |
| 6.183:DM6 | no data | - | no observed hot thin-air world fixes the box of class Z; the edges are the gap in the class table, filled and checked by the shadow sweep |
| 6.183:DM7 | no data | - | no observed methane-ethane or sulfuric-acid world is in hand beyond Titan and Venus; the ranges of U and Y are read from those two and widened |
| 6.183:DM8 | no data | - | no ammonia share is measured on an icy body outside the Solar System; the 0.03 shoulder around the 0.1 line of Class 8 is a width of the fit |
| 6.183:DM9 | no data | - | no ammonia-water ocean is observed; the box of class X is the ideal band of the ammonia-water phase diagram, read from the lab freezing curve |
| 6.91:SS16 | no data | - | the Blue edge and the Blue share are a lenient cut set on a sweep of generated bodies (Boss: near 0.8 livable); a survey of what people call a comfortable world would replace them |
| 6.91:SS17 | no data | - | the credit for a place that needs a suit is a design choice so a suit-everywhere world reads Yellow, not Red; Boss's decision on that reading would replace it |
| 6.170:RNG2 | no data | - | no measured slope of a dust ring (Jupiter, Uranus, Neptune) is in hand; a dark red-brown ring was wanted and the slope is a picture choice |
| 6.170:RNG5 | no data | - | no measured albedo of a dust ring is in hand; an educated guess between the ice ring and the dark sooty ones |
| 6.170:RNG6 | no data | - | display choices (gain, opacity floor and ceiling) of the ring picture; a rendering model that works the light of the ring would replace them |
| 6.170:RNG8 | no data | - | the drawn gap widths and the ringlet pattern are a picture; a ring structure worked out from the ring's dynamics would replace them |
| 6.170:RNG10 | no data | - | the glow of a hot body is tuned by eye (1000 K dull red, 2000 K bright orange); a rendering model of the emitted light would replace it |
| 6.170:RNG11 | no data | - | the band count is fitted to the 30 and 17 jets seen on Jupiter and Saturn; the report's jet-speed law has no usable jet speed, so no derived count is in hand |
| 6.170:RNG13 | no data | - | no source gives a jet count for Uranus or Neptune; the 5 or 7 bands are needed output |
| 6.170:RNG18 | no data | - | no source gives a storm count law by rotation or heat flux; the 1 to 4 large storms are the design's |
| 6.170:RNG19 | no data | - | band contrast and colour are tuned by eye (a methane giant nearly featureless, an ammonia giant strongly banded); the optics of the generated clouds (GEN.422 k-tables) would replace them |
| `6.3a:`--rogues`` | no data | - | the number behind each word (none, few, normal, many) is a lever of the user; the observed 2 to 10 rogues per star bounds it, and a measured rogue count would replace the ceilings |
| 6.174:SKB3 | no data | - | the ends of the faint-light table are invented to close it; the share worked from the luminosity function of the generator's own stars would replace them |
| 6.174:SKB4 | no data | - | the cell, cube map and ray sizes are a picture's needs, chosen so the band is smooth; no source gives them |
| 6.174:GLW1 | no data | - | the glow's colour is a blackbody of 4500 K because the band is worked as one V surface brightness and the page needs three channels; the mean colour of the unresolved stars by direction, from the faint-light table's stellar mix, would replace it |
| 6.176:SKT1 | no data | - | the visit budget is a serving limit chosen so a dust-free sky drops nothing; a coarser block draw for the sparse luminosity bands (a new draw and an epoch, not yet filed) would remove the need for it |
| 6.176:SKT4 | no data | - | the near-ring mass edges are spaced by hand over the lightest stars the sampler draws; a spacing taken from the mass-luminosity relation would replace them |
| 6.176:SKT5 | no data | - | the ring yield is a cost set by hand (100 stars drawn for each one that shows, as the heavier rings do); a cost per visible star set by the queue would replace it |
| 6.176:SKT6 | no data | - | the local mix weighs the young, intermediate and old populations alike, as near the Sun; the populations' densities at the observer, which a caller can pass, would replace it |
| 6.180:SKD1 | no data | - | the soft knee of the glow image is chosen by eye so the band shows at the faint end and the core does not clip; a display model of the eye's response to surface brightness would replace it |
| 6.180:SKD2 | no data | - | the grain of the glow is a seeded noise sized by eye; the Poisson variation of the stars under the limit in a cell would replace it |
| 6.180:SKD4 | no data | - | the sprite glow ramp is chosen by eye, like the stars' own scale; the same display model as SKD1 would replace it |
| 6.180:SKD5 | no data | - | the stored model ties no orbital plane to the galactic frame, so the Milky Way's plane is drawn in the plane of the orbits; a stored orientation of each system would replace it |
| 6.181:SKL1 | no data | - | the peak magnitude of the star look is chosen by eye so the faint end shows and the bright end does not clip; a display model of the eye's response fitted to the dark-adapted sky would replace it |
| 6.181:SKL2 | no data | - | the width of a drawn star and its growth per magnitude are chosen by eye so a bright star looks brighter on a monitor that cannot make it so; the eye's point-spread function would replace them |
| 6.181:SKL3 | no data | - | the colour a faint star keeps is chosen by eye; the eye's loss of colour vision at low light (the display model of SKL1) would replace it |
| 6.181:SKL4 | no data | - | the phone's limiting magnitude is set from the plan's count of about 9,000 stars in step with the 15,000 cap; the frame time a phone measures would replace it |
| 6.181:SKL5 | no data | - | the star caps are the plan's phone count and a laptop's frame time at 60 frames a second; the count the frame time of the page's own machine allows would replace them |
| 6.181:SKL6 | no data | - | the speed that takes a third light-time pass is where the residual of two passes reaches a pixel of the finest view; no measured star needs a lower one |
| 6.181:SKL7 | no data | - | the pick reach is chosen so a pixel-sized star can be hit by a mouse; a reach set by the pointer's kind (mouse, pen, touch) would replace it |
| 6.181:SKL8 | no data | - | the cap on charted names bounds the size of the route's meta; a name carried in the packed list or looked up on a pick would remove the need for it |
| 6.181:SKL9 | no data | - | the poll interval, poll count and phone width are serving limits so a poll does not hammer the route and ends; the build's own notice of completion would replace the poll |
| 6.178:PNE5 | no data | - | the hot-bottom-burning carbon to nitrogen share is set so a 5 solar mass progenitor is Type I; no source gives the surface abundances of hot-bottom-burning models by mass and metallicity |
| 6.178:PNE7 | no data | - | a safety limit that ends the carbon-rich draw loop; no draw is refused at solar metallicity |
| 6.179:NBC6 | no data | - | the grain kinds and their shares are recalled from Milky Way dust models (Draine and Li, Zubko et al.) and fitted inside the metals; no measured dust composition by nebula class is in hand |
| 6.179:NBC9 | no data | - | the dust survival shares are recalled from sputtering and reverse-shock literature; no measured dust share of an old remnant is in hand |
| 6.179:NBC10 | no data | - | the ionized, atomic and molecular shares by class are recalled from observed phase fractions; none is measured per class here |
| 6.179:NBC13 | no data | - | the wind share of a Wolf-Rayet ring is argued from the observed rings, not measured |
| 6.190:CW2 | no data | - | no census of another system's giant scatterings is in hand, so one thrown giant pair stands for the instability the hot share was set for; it goes when the scatter step takes the events' number and masses |
| 6.192:SEP1 | no data | - | L3, hand-tuned against a Hill-sphere packing sweep (125, 27 and 27 sectors at 8,128, 4,000 and 1,500 pc); it goes when a nearest-neighbour distribution from catalogued stars is fitted |
| 6.189:HC2 | no data | - | L4, a body that fits no class takes the nearest of the ordinary families by the geometric mean of its memberships with a floor of 0.01; no survey says what such a body is called, so the floor and the family list are fitted to keep the stored class defined |
| 6.18:LF63 | no data | - | the edges of the Temperature score are chosen (the Blue cut by Boss, the others proposed) and no rated catalogue of worlds exists to fit them |
| 6.18:LF64 | no data | - | no rated catalogue of worlds exists to fit the least Blue share |
| 6.18:LF65 | no data | - | no rated catalogue of worlds exists to fit the credit of a Yellow sample |
| 6.193:AW1 | no data | - | L3, the class F gravity edge 1.4 g was chosen from a sweep of 1,448 bodies; it goes when the gravity edges of classes Z and V are fixed from heavy-world physics (GEN.1136) |
