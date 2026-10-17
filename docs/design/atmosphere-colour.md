# Atmosphere colour: one consistent atmosphere, then its colour

**Status (2026-10-16):** plan; nothing is built. Boss (11:24 PM EDT, 11:25 PM EDT) made this the seminal feature of 9.0: the colour of the air of any world, worked out from the atmosphere it really has. GEN.408 checks the sources and GEN.409 turns this note into the checked design. The base is [sky-colour-analysis.md](sky-colour-analysis.md).

## 1. What Boss decided

- Store the albedo and the greenhouse values used in generation; the migration works the formula back from the stored temperature and pressure, or draws from the planet class where that fails.
- Break the atmosphere's molar density into per-gas partial pressures and calculate the molar density from them, so both match and the light meets a known mixture.
- Build ozone into the system from real research; if the research is missing, build a generator from the oxygen and the star's UV; Boss researches on request.
- Sky colour and atmosphere colour are the same; the wording is "atmosphere colour"; it replaces the planet colour where that is what it meant.
- The colour is absolute, not adapted to the star's light (card, earlier that evening), and gas giants are in the first version.
- Haze, dust and similar scatterers get reasonable draws or derivations.
- One table of every gas and compound that can occur, including every compound the descriptions mention, so that atmosphere text, composition, density and colour agree: "all of the values consistent across the board".
- Nothing is fixed per class and nothing is copied from the Solar System; every input is a stored physical quantity or a seeded draw.

## 2. The chain

| Where | Item | What |
|---|---|---|
| 8.4 phase 1 | GEN.408 | Verify the sources and say what is missing |
| 8.4 phase 1 | GEN.409 | The checked design (this note) |
| 8.4 phase 2 | GEN.410 | The gas and compound table |
| 8.4 phase 2 | GEN.411 | Molar density from the partial pressures |
| 8.4 phase 2 | GEN.412 | Albedo and greenhouse multiplier recorded, with an inverse |
| 8.4 phase 2 | GEN.422 | Stop the recalculation redraw; order the steps so pressure, mixture and temperature agree after the freeze-out |
| 8.4 phase 2 | GEN.423 | Atmosphere density, scale height and the greenhouse multiplier follow from the mixture |
| 8.4 phase 2 | GEN.424 | The balance gas is a named gas with its own molar mass; habitability reads the stored mean molar mass |
| 8.4 phase 2 | DB.67 | Columns and migration for albedo, greenhouse multiplier, molar density |
| 8.4 phase 3 | GEN.413 | Ozone from the oxygen, the star's UV and the Chapman balance |
| 8.4 phase 3 | GEN.414 | Haze, dust and cloud optical depths |
| 8.4 phase 4 | DB.68 | Columns for ozone and the aerosol depths |
| 8.4 phase 4 | GEN.415 | Description text written from the stored mixture |
| 8.4 phase 4 | TEST.209 | One consistency sweep over every class |
| 8.4 phase 4 | DOC.63 | The data model docs |
| 9.0 phase 1 | GEN.416 | Blackbody colour of a star (first part of the photometry library) |
| 9.0 phase 1 | GEN.417 | Rayleigh optical depth of the mixture |
| 9.0 phase 1 | GEN.418 | Atmosphere colour of a solid world |
| 9.0 phase 1 | GEN.419 | Atmosphere colour of gas and ice giants at a reference level |
| 9.0 phase 2 | GEN.420 | The description row and swatch in every output |
| 9.0 phase 2 | MAP.244 | The maps use the atmosphere colour in place of the planet colour |
| 9.0 phase 2 | DOC.64 | Docs of the colour |
| 9.0 phase 2 | GEN.421 | The headline item |

The foundations sit in 8.4 (four phases, each blocking the next) because the planet classes are reworked there (GEN.90) and the new classes should be born on the one table and the one mixture; the migrations take the revisions after the 0102 to 0105 stack of 8.1. VIEW.11 extends GEN.416 and stays in 9.1 with the view from a planet.

The read-only trace of 2026-10-17 (`/mnt/project-files/reports/composition-density-link-2026-10-17.md`) found that only the pressure flows from the stored mixture's side: the drawn molar density, atmosphere density, scale height and greenhouse multiplier never see the mix (molar density off by more than 10 percent for classes A, B, E, G, H, L and Q), the recalculation after orbit validation redraws the albedo, greenhouse multiplier, redox and mix, the freeze-out lowers the pressure after the temperature is set, and habitability uses 28 g/mol for every balance gas. GEN.411, GEN.422, GEN.423 and GEN.424 carry those fixes.

## 3. Open points for the design

1. **Giants:** the reference level for a gas or ice giant (cloud top or a fixed bar) and how the cloud depth sets it.
2. **Thick air:** the multiple-scattering step above tau of about 0.5 (a two-stream pass or the saturated limit).
3. **Inverse of the temperature formula:** whether the albedo and the greenhouse multiplier can both be worked back from the stored temperature and pressure (the temperature is one relation between them; whether the pressure fixes the multiplier is for GEN.412 to check in the code). Where only one can be solved, the migration draws the other from the class range, seeded.
4. **Ozone:** Earth's column near 300 Dobson units is the test target, not a preset; the oxygen and UV set it.
5. **Aerosols:** dust has no data today; its rule is a flagged parametrisation with a physical range.
6. **Stars:** the blackbody is the light; M dwarfs and cooler get a model-atmosphere table later.

## 4. Research needed

Verified: Betremieux and Kaltenegger 2014 (arXiv:1507.02107) for the refractivities of N2, CO2, H2 and He. To check (GEN.408): Sneep and Ubachs 2005, Bates 1984, Peck and Huang 1977, Penndorf 1957 (Rayleigh); Bucholtz 1995 (the Earth test value); Burkholder and Talukdar 1994, Serdyuchenko 2014, Kasting and Donahue 1980, Segura 2003 (ozone); Lemmon 2004, Tomasko 2008, Karkoschka 1994 (dust, haze, giants); Sudarsky 2000 (giant cloud classes). The study of 2026-10-17 lists how well each was checked.
