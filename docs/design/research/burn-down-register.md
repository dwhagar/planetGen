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
| 6.13:TH5 | replace | GEN.842 | Io is the one measured body; the Maxwell rheology of the tides chain (GEN.954 and GEN.955) replaces both fitted constants |
| 6.13:TH7 | replace | GEN.842 | Io is the one measured body; the Maxwell rheology of the tides chain (GEN.954 and GEN.955) replaces both fitted constants |
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
| 6.23:TR10 | no data | - | the quality factor of cold rock is measured only for the Earth and the Moon (Dziewonski and Anderson 1981, Nakamura 1982, recalled); the ceiling keeps the planets' Q inside that |
| 6.26:GS4 | replace | GEN.979 | the distribution of the last impactor by planet mass from Kokubo and Ida 2007 or Agnor et al. 1999 replaces the isolation mass and steepness; the planet-formation train can give the isolation mass of the disc (formation.isolation_mass_kg) |
| 6.26:GS5 | replace | GEN.979 | the same N-body distribution replaces the log-uniform range of the last giant impact |
| 6.26:GS6 | replace | GEN.979 | the same N-body distribution (the planetesimal-built embryos) replaces the log-uniform range |
| 6.36:PA1 | no data | - | no model of the outer disc solids is filed; the numbers are tuned to the giants at 2 to 8 AU (Fulton et al. 2021), and a measured pebble-flux model replaces them |
| 6.36:PA2 | no data | - | the share of cores that find gas is tuned to the giant occurrence and its rise with [Fe/H] (Johnson et al. 2010); no survey measures the disc lifetime it stands for |
| 6.36:PA3 | no data | - | the starting count of the inner chain is tuned to the Zhu et al. 2018 multiplicities after the formation-era dynamics build (GEN.730); the starting count of a real disc is not observed |
| 6.36:PA4 | no data | - | the masses of the Kepler planets are measured for few; the scale is tuned to the 1 to 10 Earth-mass range of Zhu and Dong 2021 and its trend with host mass and metallicity |
| 6.36:PA5 | no data | - | Weiss et al. 2018 measure the pair ratio, not the spread of masses a system starts with; tuned to the 65 percent outer-larger share |
| 6.36:PA6 | no data | - | the starting spacing is tuned to the Pu and Wu 2015 and Weiss et al. 2018 spacings after the instability build (GEN.730); the spacing before the instability is not observed |
| 6.36:PA7 | no data | - | the inner edge of the chain is measured as a fall of occurrence inside 10 days (Mulders et al. 2015); the median period and width are tuned to the Kepler rise with period |
| 6.36:PA9 | no data | - | the hot Jupiter share of giants is tuned to 0.6 to 1.2 percent of Sun-like stars and 0.27 percent of early M dwarfs; no survey measures the migration channel |
| 6.36:PA10 | no data | - | the inward move of other giants is tuned to the 1 to 3 AU peak of Fulton et al. 2021 and Cumming et al. 2008; no survey measures the migration |
| 6.36:PA11 | no data | - | the lightest oligarch that counts as a planet is tuned to the Neptune-mass counts beyond the snow line; the microlensing mass function below 1 Earth mass is unmeasured |
| 6.36:PA14 | no data | - | a requested planet class takes its mass log-uniform over 0.05 to 1000 Earth masses after 8 free draws; no survey measures the mass of a planet picked by class, so the range is a design choice kept wide |
| 6.40:RD1 | no data | - | no atmosphere of a runaway-greenhouse planet has been measured; TRAPPIST-1 b and c only forbid keeping all of it |
| 6.40:RD2 | no data | - | no abiotic oxygen air has been detected or ruled out on an M-dwarf habitable-zone planet; Luger and Barnes fix only the possible amount |
| 6.40:RD3 | no data | - | no measurement of rocky exoplanet pressures by depth below the shoreline exists; GEN.483 derives the caps at these edges from the Solar System |
| 6.40:RD5 | no data | - | the 100 K is the width over which the equilibrium lines move with pressure and metallicity; a population of hot-giant chemistry regimes has not been measured |
| 6.41:TK2 | no data | - | the line between hydrated and primitive bodies is fitted to the C to P change of the belt; no model of aqueous alteration against ice content reaches it, and no belt outside the Solar System is resolved into taxa |
| 6.41:TK5 | no data | - | the accretion spread comes from the chondrite and iron meteorite parent body ages recalled from the Solar System; no other system dates its planetesimals |
| 6.41:TK6 | replace | GEN.475 | the within-family shares are fixed fractions; the meteorite-group map of GEN.475 gives the groups reachable from the bulk composition and replaces them |
| 6.43:BM2 | no data | - | the words dense, sparse and typical are labels for the spread of belt masses at one age; no survey resolves the masses of belts of asteroids, and the factor of 3 is a choice of where the words change |
| 6.44:MG2 | no data | - | the widths of the axes and the distance past which a body is called none are choices; the collection's spread inside each group is not in the table, which gives averages and ranges only |
| 6.47:SE4 | no data | - | the tail slope of the largest solar particle events is fitted to one Carrington-class event in 500 years; a longer proxy record of the largest events (nitrates, cosmogenic isotopes) would replace it |
| 6.14:SP11 | replace | GEN.979 | Kokubo and Genda 2010 equation 4 with the distribution of the last impacts' mass ratio (GEN.978, `impact_spin`) replaces the interim scale |
| 6.14:SP12 | replace | GEN.979 | the chance of a giant impact by planet mass (`impact_spin.giant_impact_probability`, GEN.978) replaces the interim small-planet factor |
| 6.51:PS14 | no data | - | the label is a margin between the model's own solar-metallicity and metal-poor stars; no survey count of subdwarfs by luminosity deficit is held, and the first one replaces the margin |
| 6.60:GM2 | replace | GEN.1006 | Mustill and Villaver 2012 (read) print three Jovian points and a terrestrial bracket; the multi-star build replaces the linear lookup with an interpolation through them |
| 6.60:GM3 | no data | - | a design default of the multi-star build; no measured limit for a close pair's separation in summed radii is held, and an N-body audit (GEN.1006 section 9) is optional; tier T6 stays |
| 6.60:GM7 | replace | GEN.1006 | the 3+1 network of Vynatheya, Mardling and Hamers 2023 replaces the margin inside its training range; outside it the margin stays, flagged lower confidence |
| 6.60:GM12 | replace | GEN.1028 | the criterion holds only for q_out up to 5; GEN.1028 enforces the range and GEN.1016 replaces it with the quadruple fit and the 3+1 criterion beyond it |
| 5.2:BC20 | replace | GEN.1057 | the conditional ring draw (a ring only after a moon fell or a comet passed inside the Roche limit, `ring_origin`) replaces the kind weights and the chance of a ring; the Solar System's four giants are one system, so the weights have no second sample |
| 5.2:BC21 | no data | - | the edges of ring systems are measured for four giants only; the outer edge against the Roche limit and the inner edge in planet radii have no second sample |
| 5.2:BC22 | replace | GEN.1057 | the ring's mass follows from the disrupted moon's or comet's mass and its spreading (`ring_origin`) once the conditional draw is wired |
| 5.2:BC23 | replace | GEN.1055 | the ice share and darkness follow from the ring's age and pollution by the snow line (GEN.1055) |
| 6.70:RO1 | no data | - | no simulation table of the ring left by a disrupted moon by composition is in hand; Canup 2010 and Black and Mittal 2015 disagree on the share |
| 6.70:RO2 | no data | - | the age of Saturn's rings is contested between 1e8 years and the age of the planet, and no other ring system is dated |
| 6.70:RO3 | no data | - | no survey measures how often a ring-making comet passes inside a giant's Roche limit |
| 6.70:RO5 | no data | - | no N-body table of the share of a torn comet that stays in orbit by speed and periapsis is in hand; Dones 1991 gives a range only |
