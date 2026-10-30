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
| 6.7:AM81 | replace | GEN.722 | the radius valley follow-up settles the mass-loss mechanism and the widths; Kepler occurrence tables are the data (Fulton 2017, Van Eylen 2018 as second-hand widths) |
| 6.10:RV2 | replace | GEN.722 | the radius valley follow-up settles the mass-loss mechanism and the widths; Kepler occurrence tables are the data (Fulton 2017, Van Eylen 2018 as second-hand widths) |
| 6.10:RV3 | replace | GEN.722 | the radius valley follow-up settles the mass-loss mechanism and the widths; Kepler occurrence tables are the data (Fulton 2017, Van Eylen 2018 as second-hand widths) |
| 6.9:SP1 | no data | - | a display choice (how fast a rendered body turns on screen), not a physical quantity; tier T9 stays |
| 6.13:TH5 | replace | GEN.842 | Io is the one measured body; the Maxwell rheology of the tides chain (GEN.954 and GEN.955) replaces both fitted constants |
| 6.13:TH7 | replace | GEN.842 | Io is the one measured body; the Maxwell rheology of the tides chain (GEN.954 and GEN.955) replaces both fitted constants |
| 6.15:RB6 | no data | - | AP-9, AE-9 and the Divine-Garrett model are not held and the Galileo and Juno dose records are figures; the shapes stay recalled stand-ins |
| 6.15:RB7 | no data | - | AP-9, AE-9 and the Divine-Garrett model are not held and the Galileo and Juno dose records are figures; the shapes stay recalled stand-ins |
| 6.18:LF5 | no data | - | the pace timelines are scaled by hand from one history; the hard-step life timeline (GEN.750, closed) supersedes them for life and the table goes when the civilization check stops reading it |
| 6.18:LF6 | no data | - | design thresholds of the life scores (GEN.649): no measured quantity exists; tier T9 stays |
| 6.18:LF7 | no data | - | design thresholds of the life scores (GEN.649): no measured quantity exists; tier T9 stays |
| 6.18:LF8 | no data | - | design thresholds of the life scores (GEN.649): no measured quantity exists; tier T9 stays |
| 6.18:LF9 | replace | GEN.894 | the strange-solvent classes replace the solvent factors and the methane-sea gate; no measured limit of a non-water biochemistry exists |
| 6.18:LF24 | replace | GEN.894 | the strange-solvent classes replace the solvent factors and the methane-sea gate; no measured limit of a non-water biochemistry exists |
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
| 6.19:LF40 | replace | GEN.894 | the strange-solvent classes 5 to 9 carry placeholder supplies and windows; there is no organism of these classes to measure |
| 6.19:LF41 | replace | GEN.894 | the strange-solvent classes 5 to 9 carry placeholder supplies and windows; there is no organism of these classes to measure |
| 6.19:LF42 | replace | GEN.894 | the strange-solvent classes 5 to 9 carry placeholder supplies and windows; there is no organism of these classes to measure |
| 6.19:LF43 | replace | GEN.894 | the strange-solvent classes 5 to 9 carry placeholder supplies and windows; there is no organism of these classes to measure |
| 6.19:LF44 | replace | GEN.894 | the strange-solvent classes 5 to 9 carry placeholder supplies and windows; there is no organism of these classes to measure |
| 6.19:LF45 | replace | GEN.894 | the strange-solvent classes 5 to 9 carry placeholder supplies and windows; there is no organism of these classes to measure |
| 6.19:LF47 | no data | - | chosen or Earth-calibrated numbers of the life supplies and ledgers (each carries its own tier in its file); no other biosphere exists to calibrate them |
| 6.19:LF59 | no data | - | chosen or Earth-calibrated numbers of the life supplies and ledgers (each carries its own tier in its file); no other biosphere exists to calibrate them |
| 6.19:LF60 | no data | - | chosen or Earth-calibrated numbers of the life supplies and ledgers (each carries its own tier in its file); no other biosphere exists to calibrate them |
| 6.19:LF53 | no data | - | chosen or Earth-calibrated numbers of the life supplies and ledgers (each carries its own tier in its file); no other biosphere exists to calibrate them |
| 6.20:PM1 | replace | GEN.728 | the observed occurrence rates (Fulton 2021, Petigura 2018, Dressing and Charbonneau 2015) enter as bands the planet count and spacing items (GEN.728, GEN.729, GEN.828, GEN.829) fit to |
| 6.20:PM2 | replace | GEN.728 | the observed occurrence rates (Fulton 2021, Petigura 2018, Dressing and Charbonneau 2015) enter as bands the planet count and spacing items (GEN.728, GEN.729, GEN.828, GEN.829) fit to |
| 6.20:PM6 | replace | GEN.728 | the observed occurrence rates (Fulton 2021, Petigura 2018, Dressing and Charbonneau 2015) enter as bands the planet count and spacing items (GEN.728, GEN.729, GEN.828, GEN.829) fit to |
| 6.20:PM7 | replace | GEN.728 | the observed occurrence rates (Fulton 2021, Petigura 2018, Dressing and Charbonneau 2015) enter as bands the planet count and spacing items (GEN.728, GEN.729, GEN.828, GEN.829) fit to |
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
| 6.32:SE4 | no data | - | the rate of events above 10 000 pfu rests on one recalled count of a Carrington-class event (McCracken 2001); the NOAA list ends at 9 events above 10 000 pfu in 42 years, and no proxy record of the largest events (nitrates, cosmogenic isotopes) is held or reachable; any measurement replaces the slope |
