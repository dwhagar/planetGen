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
| 6.7:AM64 | replace | GEN.804 | only figure read-offs of Johnstone et al. 2018 exist; the upper-atmosphere model that reads the stored XUV replaces the fitted exobase constants |
| 6.7:AM65 | replace | GEN.804 | only figure read-offs of Johnstone et al. 2018 exist; the upper-atmosphere model that reads the stored XUV replaces the fitted exobase constants |
| 6.7:AM66 | replace | GEN.804 | only figure read-offs of Johnstone et al. 2018 exist; the upper-atmosphere model that reads the stored XUV replaces the fitted exobase constants |
| 6.7:AM67 | replace | GEN.804 | only figure read-offs of Johnstone et al. 2018 exist; the upper-atmosphere model that reads the stored XUV replaces the fitted exobase constants |
| 6.7:AM71 | replace | GEN.804 | only figure read-offs of Johnstone et al. 2018 exist; the upper-atmosphere model that reads the stored XUV replaces the fitted exobase constants |
| 6.7:AM81 | no data | - | the radius valley follow-up (GEN.722) is retired; the widths are fitted to the Kepler occurrence tables (Fulton 2017, Van Eylen 2018 as second-hand widths) and any measured mass-loss width replaces them |
| 6.10:RV2 | no data | - | the radius valley follow-up (GEN.722) is retired; the widths are fitted to the Kepler occurrence tables (Fulton 2017, Van Eylen 2018 as second-hand widths) and any measured mass-loss width replaces them |
| 6.10:RV3 | no data | - | the radius valley follow-up (GEN.722) is retired; the widths are fitted to the Kepler occurrence tables (Fulton 2017, Van Eylen 2018 as second-hand widths) and any measured mass-loss width replaces them |
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
| 6.18:LF23 | no data | - | the animal oxygen floor by grade is a literature number (Mills 2014, recalled) for one grade; GEN.877 is closed and kept the band; the widths are chosen |
| 6.18:LF25 | no data | - | one score per rung of Table 2 of the Exoplanetary Habitability doc; design numbers with no data |
| 6.18:LF49 | no data | - | the hard-step means follow Carter 1983 and Watson 2008 in form; the only calibration is one planet's history, and the rule against fixed Solar System content keeps it a check |
| 6.18:LF50 | no data | - | the hard-step means follow Carter 1983 and Watson 2008 in form; the only calibration is one planet's history, and the rule against fixed Solar System content keeps it a check |
| 6.18:LF51 | no data | - | the hard-step means follow Carter 1983 and Watson 2008 in form; the only calibration is one planet's history, and the rule against fixed Solar System content keeps it a check |
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
| 6.19:LF61 | no data | - | the gates and the liquid window of the acid cloud are chosen limits and the acid strength line is fitted through two Venus cloud levels; no second acid world exists to calibrate against, and the Zeleznik 1991 curve is not in the reachable papers |
| 6.19:LF62 | no data | - | the ammonia share of an icy body's ices is a choice between the comet value and the cosmic cap; no survey measures the ammonia of the ocean of any world but a few plume spectra, and the Mg/Si line is the textbook one |
| 6.19:LF53 | no data | - | chosen or Earth-calibrated numbers of the life supplies and ledgers (each carries its own tier in its file); no other biosphere exists to calibrate them |
| 6.20:PM2 | no data | - | the observed occurrence rates (Fulton 2021, Petigura 2018, Dressing and Charbonneau 2015) are the bands the planet count and spacing fit to (GEN.728, GEN.729, GEN.828, GEN.829, built); the steepening itself is tuned to them and no disc survey measures it |
| 6.21:BO1 | no data | - | the oxygen balance has one measured world; Earth's fluxes are recalled (Hayes and Waldbauer 2006, Holland 2002) and calibrate the scale; tier stays |
| 6.21:BO2 | no data | - | the oxygen balance has one measured world; Earth's fluxes are recalled (Hayes and Waldbauer 2006, Holland 2002) and calibrate the scale; tier stays |
| 6.21:BO3 | no data | - | the oxygen balance has one measured world; Earth's fluxes are recalled (Hayes and Waldbauer 2006, Holland 2002) and calibrate the scale; tier stays |
| 6.21:BO4 | no data | - | the oxygen balance has one measured world; Earth's fluxes are recalled (Hayes and Waldbauer 2006, Holland 2002) and calibrate the scale; tier stays |
| 6.21:BO6 | no data | - | the oxygen balance has one measured world; Earth's fluxes are recalled (Hayes and Waldbauer 2006, Holland 2002) and calibrate the scale; tier stays |
| 6.21:BO8 | no data | - | the oxygen balance has one measured world; Earth's fluxes are recalled (Hayes and Waldbauer 2006, Holland 2002) and calibrate the scale; tier stays |
| 6.22:CC1 | replace | GEN.962 | the moon rate item (giant-impact and capture frequency) replaces the moon mass floor and function; the satellite mass ratios are a Solar System check, not a source |
| 6.22:CC2 | replace | GEN.962 | the moon rate item (giant-impact and capture frequency) replaces the moon mass floor and function; the satellite mass ratios are a Solar System check, not a source |
| 6.22:CC4 | replace | GEN.684 | albedo and greenhouse from the composition (the ground and cloud albedo) delete the regime ranges |
| 6.22:CC6 | no data | - | the lightest planet is a naming boundary and the habitable-request steering is a user-interface choice; neither is a measured quantity |
| 6.22:CC10 | no data | - | the lightest planet is a naming boundary and the habitable-request steering is a user-interface choice; neither is a measured quantity |
| 6.23:TR3 | no data | - | the form follows from physics (a suspension has no static rigidity); the floor is the lowest value that keeps the Io state's k2 under the fluid limit, and no measurement of the rigidity of a partial melt at planetary scale exists |
| 6.23:TR10 | no data | - | measured tidal Q of cold rock (Moon 38 to 41, Mars 86) lies far under the ceiling, so it constrains the ceiling only from below; the seismic Q the row cited is another frequency; no tidal Q above 100 is measured; tier T6, tweak L3 stay |
| 6.26:GS4 | no data | - | the distribution of the last impactor by planet mass is in Kokubo and Ida 2007 and Agnor et al. 1999 as figures and fits we do not hold as a table; a table of the last impact's mass ratio by planet mass replaces the isolation mass and steepness (the planet-formation train can give the isolation mass of the disc, formation.isolation_mass_kg) |
| 6.26:GS5 | no data | - | the same N-body table replaces the log-uniform range of the last giant impact |
| 6.26:GS6 | no data | - | the same N-body table (the planetesimal-built embryos) replaces the log-uniform range |
| 6.36:PA1 | no data | - | no model of the outer disc solids is filed; the numbers are tuned to the giants at 2 to 8 AU (Fulton et al. 2021), and a measured pebble-flux model replaces them |
| 6.36:PA2 | no data | - | the share of cores that find gas is tuned to the giant occurrence and its rise with [Fe/H] (Johnson et al. 2010); no survey measures the disc lifetime it stands for |
| 6.36:PA3 | no data | - | the starting count of the inner chain is tuned to the Zhu et al. 2018 multiplicities after the formation-era dynamics build (GEN.730); the starting count of a real disc is not observed |
| 6.36:PA4 | no data | - | the masses of the Kepler planets are measured for few; the scale is tuned to the 1 to 10 Earth-mass range of Zhu and Dong 2021 and its trend with host mass and metallicity |
| 6.36:PA5 | no data | - | Weiss et al. 2018 measure the pair ratio, not the spread of masses a system starts with; tuned to the 65 percent outer-larger share |
| 6.36:PA6 | no data | - | the starting spacing is tuned to the Pu and Wu 2015 and Weiss et al. 2018 spacings after the instability build (GEN.730); the spacing before the instability is not observed |
| 6.36:PA7 | no data | - | the inner edge of the chain is measured as a fall of occurrence inside 10 days (Mulders et al. 2015); the median period and width are tuned to the Kepler rise with period |
| 6.36:PA9 | no data | - | the share 0.0092 now sits between two read statistics (Zhu and Dong 2021: 0.62 percent Kepler, 0.9 to 1.2 percent radial velocity), but the row also holds the metallicity exponent 2.2 and the period range, which no survey constrains separately; the row stays L3 until they are split |
| 6.36:PA10 | no data | - | the inward move of other giants is tuned to the 1 to 3 AU peak of Fulton et al. 2021 and Cumming et al. 2008; no survey measures the migration |
| 6.36:PA11 | no data | - | Suzuki et al. 2016 (read) give the planet counts per log mass ratio beyond the snow line down to their sensitivity edge but no floor below it; the lightest oligarch that counts stays tuned to those counts; tier T5, tweak L3 stay |
| 6.36:PA14 | no data | - | a requested planet class takes its mass log-uniform over 0.05 to 1000 Earth masses after 8 free draws; no survey measures the mass of a planet picked by class, so the range is a design choice kept wide |
| 6.40:RD1 | no data | - | no atmosphere of a runaway-greenhouse planet has been measured; TRAPPIST-1 b and c only forbid keeping all of it |
| 6.40:RD2 | no data | - | no abiotic oxygen air has been detected or ruled out on an M-dwarf habitable-zone planet; Luger and Barnes fix only the possible amount |
| 6.40:RD3 | no data | - | no measurement of rocky exoplanet pressures by depth below the shoreline exists; GEN.483 derives the caps at these edges from the Solar System |
| 6.40:RD5 | no data | - | the 100 K is the width over which the equilibrium lines move with pressure and metallicity; a population of hot-giant chemistry regimes has not been measured |
| 6.41:TK2 | no data | - | the line between hydrated and primitive bodies is fitted to the C to P change of the belt; no model of aqueous alteration against ice content reaches it, and no belt outside the Solar System is resolved into taxa |
| 6.41:TK5 | no data | - | the accretion spread was recalled and is now read: Lichtenberg et al. 2021 give core formation under about 1 Myr in the first reservoir and 1.8 to 4 Myr in the second, consistent with the 3 Myr; the value remains the one Solar System spread tuned by the melting test, no other system dates its planetesimals; tier T3, tweak L3 stay |
| 6.43:BM2 | no data | - | the words dense, sparse and typical are labels for the spread of belt masses at one age; no survey resolves the masses of belts of asteroids, and the factor of 3 is a choice of where the words change |
| 6.44:MG2 | no data | - | the widths of the axes and the distance past which a body is called none are choices; the collection's spread inside each group is not in the table, which gives averages and ranges only |
| 6.47:SE4 | replaced | - | the tail now comes from the isotope record of Usoskin and Kovaltsov 2012 (0.0077 events a year at 1e10 cm-2, power law -4 beyond); no longer hand-tuned; the pairing of fluence with peak flux stays McCracken's recalled number (tier T3, tweak L1) |
| 6.14:SP22 | no data | - | the generator holds no secular frequencies or fossil shapes to give the Cassini-state share of the inclination; a Laplace-Lagrange run of the generated system with the body's gravity field replaces the log-uniform share, which is kept inside Mercury's 0.004 and the bound 1
| 6.51:PS14 | no data | - | the label is a margin between the model's own solar-metallicity and metal-poor stars; no survey count of subdwarfs by luminosity deficit is held, and the first one replaces the margin |
| 6.60:GM3 | no data | - | a design default of the multi-star build; no measured limit for a close pair's separation in summed radii is held, and an N-body audit (GEN.1006 section 9) is optional; tier T6 stays |
| 6.60:GM7 | replace | GEN.1006 | the 3+1 network of Vynatheya, Mardling and Hamers 2023 replaces the margin inside its training range; outside it the margin stays, flagged lower confidence |
| 6.60:GM12 | replace | GEN.1016 | the criterion holds only for q_out up to 5; GEN.1028 enforces the range (built), and GEN.1016 replaces it with the quadruple fit and the 3+1 criterion beyond it |
| 6.60:GM13 | replace | GEN.1006 | the 3+1 axis margin of 1.3 is hand-tuned on 331 N-body systems; a larger own N-body sample (GEN.1006 section 9) or a retrained network replaces it |
| 6.60:GM14 | replaced | - | GEN.1155: the 2+2 axis margin is now 1.5, tuned (L2) to the paper's 0.94 precision on 302 + 185 own REBOUND 2+2 systems; it was the 3+1 margin of 1.3 with no 2+2 sample |
| 6.148:RGD1 | no data | - | the clean-ice reflectance and the half-darkening share of an ice and dark mixture are recalled (Clark and Lucey 1984 was not read); the share is tuned to Saturn's ring reflectances, and a measured mixture curve replaces both |
| 6.71:RR1 | no data | - | the rate of a resonance ring on a small body has no survey behind it: the three known rings are of bodies chosen for occultation campaigns, and the count of bodies searched is not published as a rate |
| 6.71:RR4 | no data | - | the mass and width range of a thin resonance ring is bracketed by the optical depth it must give; the measured optical depth and width of Quaoar's ring were not read |
| 6.69:LQ3 | no data | - | the brine and ammonia-water eutectic densities, viscosities and boiling points are recalled values; laboratory brine tables (Pitzer-model compilations) can replace them when read |
| 6.69:LQ5 | no data | - | the molten iron sulfide and Fe-S eutectic melt properties are recalled; measured Fe-S melt properties (Terasaki et al.) replace them when read |
| 5.2:BC20 | no data | - | the share of giants with a ring is now an output of the fall and passage histories (`ring_origin`); no ring census of exoplanets is held, and Cassini's ring age and mass are argued, so the history's own stand-ins RO1 to RO6 carry it |
| 5.2:BC21 | no data | - | the edges of ring systems are measured for four giants only; the outer edge against the Roche limit and the inner edge in planet radii have no second sample |
| 5.2:BC22 | no data | - | the kind limits are the middles of the gaps between four Solar System rings; a second sample of measured ring masses replaces them |
| 5.2:BC23 | replace | GEN.1055 | the ice share and darkness follow from the ring's age and pollution by the snow line (GEN.1055) |
| 6.70:RO1 | no data | - | no simulation table of the ring left by a disrupted moon by composition is in hand; Canup 2010 and Black and Mittal 2015 disagree on the share |
| 6.70:RO2 | no data | - | the age of Saturn's rings is contested between 1e8 years and the age of the planet, and no other ring system is dated |
| 6.70:RO3 | no data | - | no survey measures how often a ring-making comet passes inside a giant's Roche limit |
| 6.70:RO5 | no data | - | no N-body table of the share of a torn comet that stays in orbit by speed and periapsis is in hand; Dones 1991 gives a range only |
| 6.74:SH3 | no data | - | the pressure a given shoreline ratio allows is not measured as a table: only Mars's 6 mbar and Pluto's 10 micro-bar sit at ratios of 1.8 and 0.78; the caps stay the design note's steps until a survey of thin atmospheres fits them |
| 6.78:GT3 | replace | GEN.1058 | the share of ice giants with a recorded merger and the impactor mass share are tuned to one Solar System count (Uranus); the impact history of the outer chain, once recorded, gives the merger and its impactor |
| 6.96:IO2 | no data | - | fitted with the Jupiter row to the retrograde to prograde split; a capture simulation's width replaces it when one is held |
| 6.96:IO4 | no data | - | no family velocity field is held beyond the recalled tens to hundreds of m/s; the first one replaces the range |
| 6.83:SX2 | no data | - | no simulation or survey gives the captured share of a swept span for a given migration; replaced by the resonant populations of a Neptune-migration run |
| 6.83:SX3 | no data | - | the ejected share of a scattering giant is not tabulated by Safronov number; replaced when the scattering simulations are read |
| 6.83:SX4 | no data | - | no survey statistic gives a hot to cold share for an arbitrary system; the observed classical populations are the only count |
| 6.83:SX5 | no data | - | the split of the bound swept mass among scattered, centaur, detached, trojan, implanted and moon populations is not measured; replaced by a population synthesis |
| 6.95:BR20 | no data | - | the size-frequency of the main belt above 40 km (Bottke et al. 2005 and the catalogue counts) was not reachable: the paper is not on arXiv and no catalogue was fetched; Ivezic et al. measured only to 40 km. Any census of bodies above 100 km replaces the break and the steep slope |
| 6.100:CM14 | no data | - | the share a gas disc adds by drag is in no paper reached: Pollack et al. 1979 is the abstract only (NTRS), Cuk and Burns 2004 is not reachable, and the two reviews read (Jewitt and Haghighipour 2007, Sheppard 2006) give the gate and size window (CM18) and call the model weakly constrained; the first measured capture simulation of the drag channel replaces it |
| 6.98:BW2 | no data | - | no measurement of surface winds across the thin-air range (Mars, Pluto, Triton) gives the fade's shape; the first drag-law fit from surface wind data replaces it |
| 6.98:BW4 | no data | - | tuned to the Earth's mid-latitude surface wind; a boundary-layer drag law on the surface roughness replaces it when the surface roughness is stored |
| 6.98:BW5 | no data | - | the Earth's night to day wind ratio is the only value held; a stable boundary-layer model replaces it |
| 6.88:RG2 | replace | GEN.1111 | the ejection records (RP1) give the time of each loss after formation, so the drawn delay distribution is replaced by the architecture's own |
| 6.87:RK7 | replace | GEN.1111 | the closest approach of each scattering is in the architecture's own ejection record, which replaces the line through three medians |
| 6.112:MO2 | replace | GEN.1111 | the ejection record gives the real closest approach, so the curve through two printed points is not needed |
| 6.86:TU1 | replace | GEN.1111 | the rogue research settles the per-star rate of free-floating planets and replaces the level ceilings |
| 6.86:TU2 | no data | - | uniform placement is the neutral choice; the smoothed rogue field (GEN.1115) is the only thing that could change it |
| 6.86:TU3 | replace | GEN.1111 | the ejection records give the scatterer's distance, so the drawn range is replaced by the architecture's own |
| 6.108:TF1 | replace | GEN.1103 | a physical source for the mass the Rogues level cannot hold (the unbound remnants of the channels of section 6.103) replaces the stand-in origin |
| 6.108:TF2 | replace | GEN.1105 | the size-frequency law of the collisional fields sets the largest body and the density of the pieces |
| 6.108:TF3 | no data | - | nothing measures the spread of a field that was never observed; the dispersion of a real source (IF12) is the only thing that could replace it |
| 6.108:TF4 | no data | - | the counts and masses are chosen to keep a sector's top-up to a dozen clouds; no measurement exists to replace them |
| 6.114:RY1 | no data | - | no stagnant-lid eruption model yet; the one-box mantle temperature does not discriminate between bodies, so the cut stays a heat flow between Mars and the Earth |
| 6.117:BHS3 | no data | - | no source gives the share of a belt (as against a planet system) that a star's sudden mass loss destabilises; the value is set so the white dwarf pollution the model gives is in the observed range, and any N-body run of a belt round a mass-losing star replaces it |
| 6.118:IP5 | no data | - | no paper gives the gas-drag share against the three-body channel (Pollack et al. 1979 only as an abstract, Cuk and Burns 2004 not reached); a gas-drag capture run of the generated envelope replaces it |
| 6.118:IP8 | no data | - | no binary-exchange run gives the size of the captured body for a general disc; the design's few hundred km is the only guide |
| 6.118:IP10 | no data | - | the disruption size ratio is an energy argument from recalled strengths, not a measurement; the size-dependent strength of a cascade run of the generator's own population replaces it |
| 6.119:CP3 | no data | - | no source gives how many pieces a captured parent was made of; the Himalia group's four known members are the only guide, and a parent built from the disc's own bodies (belt design) replaces it |
| 6.121:TJ4 | no data | - | no source read gives the leading share of a swarm; the capture history of the generated giant replaces it |
| 6.121:TJ5 | no data | - | no source read gives the eccentricity and inclination limits; an N-body map of the trojan region by mass ratio replaces them |
| 6.41:TK6 | no data | - | the within-family shares are fixed fractions; the meteorite-group map of GEN.475 labels the group of a mix but gives no shares for the taxa of a family, and no table of asteroid taxa by mass in a belt exists to fit them to |
| 6.130:AB1 | replace | GEN.804 | the exponents are tuned to the earlier table's level; the rotation and X-ray evolution of Johnstone et al. 2021 (the model, not its figure) gives the decay itself, and the upper-atmosphere model that reads the stored XUV is the item that needs it |
| 6.43:LB1 | replace | GEN.1000 | the size law of GEN.566 gives each belt's largest body; until it is built the size is a range drawn per belt |
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
| 6.129:IT6 | no data | - | the 1993 ammonia-water constants of Durham, Kirby and Stern are not in hand; they replace the three interpolated points |
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
| 6.163:IHA1 | no data | - | no measured viscosity of any rocky planet's mantle other than the Earth's is in hand; the prefactor is retuned to hold the reference viscosity, not fitted to a planet |
