# Units, number formatting and the shared unit-ladder module

How every measurement planetGen shows should be turned into text: one ladder
of units per quantity, one number rule, one place where the data lives, and
one generic algorithm in Python and in the browser. This note holds the
proposed ladders for the 14 families the repository shows today, the
module API, the migration order, the defects found in the current formatter
and the tests that keep the Python and JavaScript copies identical.

Informs: UX.23, UX.22, UX.36, UX.30, UX.32, GEN.66 (and, only where relevant, UX.3 and UX.42)

Status: research, 2026-10-09; decisions marked "Boss" are his, everything else is a recommendation

Evidence tags: [S] seen in a search result or a fetched primary file (URL in
Sources), [C] computed or measured in the research sandbox, [R] recalled and
unconfirmed. The research environment could only read search-result text and
raw files, not whole papers or standards, and the web-search budget ran out
partway through, so the screen-reader, exoplanet-catalogue and unit-toggle
points are [R]. Every [R] is on the Evidence notes list.

## Decisions already taken

- Distances (Boss, 2026-09-30): a parsec-family value carries one familiar unit in parentheses, "4.2 pc (13.7 ly)", "2.4 mpc (495 AU)" (`DISTANCE_PAREN_MIN_LY`, `DISTANCE_PAREN_MIN_AU` in `util/format.py`).
- Body radii (Boss, 2026-09-30): a planet, moon or star radius is always km in scientific notation, the one exception to the distance ladder (`format_body_radius_km`).
- Speed (Boss, 2026-10-01, UX.13): "km/h on the low speed end to Mm/s on the high speed end", then multiples of c from a tenth of light speed.
- UX.22 (Boss, 2026-10-01 15:10Z): "standardize ALL measurements into trees like we have so that we always have meaningful units ... just everything that can have units. Atmospheric pressure and surface conditions should show customary units as well as a secondary to help contextualize the metric values given." Surface temperature in K, °C and °F and pressure in Pa-ladder, atm and psi shipped in PR #234.
- Number rule (UX.20, UX.36, Boss 2026-10-02): scientific notation for a whole number only from the 7th digit, and from the 5th digit before the point when decimals are shown. Implemented and verified (section 3).
- UX.23 plan (2026-10-07): the ladder is built on astropy.units (GEN.66, done in PR #767).

## Summary and recommendation

1. **Build the module data-driven, not as more hand-written twins.** One Python table of ladders (rungs with a symbol, an SI factor, an optional start threshold, a singular form, a spelled-out name and an optional aside unit), one generic picking-and-formatting function, and a generated JS data file consumed by one generic JS function. A prototype of this reproduced the existing `format_distance_m`, `format_speed_kms` and `format_duration_seconds` exactly on 600,000 random values with 0 mismatches [C], so step 0 of the migration can change no output at all. The generated file also removes the three hand-kept JS twins (`distance.js`, `speed.js`, `period.js`, 189 lines) and the constants hard-coded in them.
2. **Use astropy.units for the numbers and the tests, not in the display path.** Derive every rung factor from `astropy.units` and `astropy.constants` at import (as `physics/constants.py` already does since GEN.66) and cross-check every rung in tests. Format plain floats: a plain divide costs 0.17 to 0.5 µs, `Unit.to(value)` 4.5 to 5.6 µs, a Quantity round trip 22 to 37 µs, pint 87 µs [C]. Speed is not decisive; astropy's output style ("km s⁻¹", no grouping, no UX.36 rule), numpy scalars and version drift are. Do not adopt pint (its 0.24.1 fails to import on Python 3.9) or a JS unit library (none knows AU, pc or Earth and Jupiter masses). Detail in section 2.
3. **Pin astropy's constant sets: a real reproducibility gap in the lock file today.** `requirements.lock` resolves astropy 6.0.1 (Python 3.9), 6.1.7 (3.10) and 8.0.1 (3.11 and later), and astropy 8.0 switched its default physical constants from CODATA 2018 to CODATA 2022 [C]. G, M_sun, M_earth, M_jup, R_sun, L_sun and the au are identical, but the proton and electron masses differ by about 1.4e-9, so `HYDROGEN_ATOM_MASS_KG` differs between a 3.9 and a 3.12 install, and `tests/test_astropy_constants.py` allows 1e-4 on it. Fix: call `astropy.physical_constants.set("codata2018")` and `astropy.astronomical_constants.set("iau2015")` once, before `astropy.units` is first imported (it raises `RuntimeError: astropy.units is already imported` otherwise [C]); this gives identical values on 6.0.1 and 8.0.1 [C]. The IAU 2015 B3 nominal values (R_sun 6.957e8 m, L_sun 3.828e26 W, R_jup 7.1492e7 m, GM_sun 1.3271244e20 m³/s²) are the same in every astropy version [C].
4. **Fix four concrete defects in the current formatter while moving it.**
   - (a) Python and JS round half-way values differently: `_three_figures(9.995)` gives "9.99" in Python and "10" in JS, and `1.005` gives "1" versus "1.01" [C]. Round half-up on the shortest decimal representation in both (`Decimal(repr(x))` with `ROUND_HALF_UP`); that matched `Intl.NumberFormat` on 400,000 cases including ties, at 5.8 µs per call against 1.3 µs for `format()` [C].
   - (b) Between 1,000,000 km and 1 AU every distance is scientific km ("5.79 × 10⁷ km" for Mercury's 0.387 AU) [C]. Start AU at 1,000,000 km (0.0067 AU).
   - (c) A value below the bottom rung prints a row of zeros ("0.00000000000025 Pa", a Hawking luminosity) [C]. Use scientific text when the mantissa is under 1e-3.
   - (d) "Gy" means gigayear in the period ladder but is the gray in SI, and GEN.87 will add doses. Use "kyr", "Myr" and "Gyr".
5. **Every quantity gets one ladder and one aside rule** (section 1). The rogue-planet and phenomenon pages (`web/system_pages.py`) and the 3D views (`systemview3d.js`, `galaxysystem.js`, `galaxymap3d.js`) still bypass the ladders; section 5 names the files.
6. **Wrap values in `<span class="qty" data-q="length" data-si="…">` and keep NBSP out of the text.** That gives a no-wrap rule, tabular numerals, a hook for a unit preference and for copy behaviour, without putting non-breaking spaces into the Markdown and wikitext the same formatter writes.

## 1. What the project displays today, and the proposed ladders

### Where the formatter lives now

Function names are given instead of line numbers, which drift; grep for them.

| Piece | File | Notes |
|---|---|---|
| Number rule (UX.20, UX.36) | `src/planetgen/util/format.py`: `format_number`, `scientific_text`, `_shows_too_many_digits`, `_three_figures` | 7+ whole digits with no decimals, or 5+ digits before the point with decimals, becomes "1.23 × 10⁶" (3 significant figures, Unicode superscripts). |
| Ladders | `format.py`: `DISTANCE_LADDER`, `SPEED_LADDER`, `PERIOD_LADDER`, `PRESSURE_LADDER`; `format_temperature_k` | Hand-written per quantity, each with its own loop. |
| HTML wrappers | `src/planetgen/web/lib/fmt.py` | `dash_unless_finite(...)` around each formatter (TEST.53). |
| Ad hoc | `web/lib/tabledisplay.py`, `web/system_pages.py` (`FIELD_SPECS`, `_temperature_text`, `_pressure_text`), `web/lib/systempage.py` (`_gravity_text`), `generation/planet.py` (gravity line) | Mass as "kg (x Earth masses)", luminosity as "W (x Sol)", Gy, G, W/m², particles/cm³ all formatted inline. |
| JS twins | `src/html/static/numberformat.js`, `distance.js`, `speed.js`, `period.js` | Tested against Python by spawning node (`tests/test_number_format.py`, `test_distance_format.py`, `test_speed_period_format.py`). |
| JS bypassing the ladders | `systemview3d.js` (orbit and period fields), `galaxysystem.js` (period), `galaxymap3d.js` (17 `toFixed`/`toPrecision` calls), `galaxystages.js` (sector counts) | `toPrecision(3)` gives "1.50e+3" style output for four-digit values; `Math.round(radius).toLocaleString()` ignores the "radius is always scientific km" rule. |

Stored strings are not a migration worry: schema v5 dropped the pre-rendered `table_*` columns, so display text is computed at read time. The Markdown and wikitext render in `generation/*.py` still bakes strings (for example `star.py` writes `f"{format_number(self.temperature)} K"`), so those files need the same wrappers.

### Existing ladder spans [C]

The existing rungs cover: distance km to 1.5e8 km, AU 1 to 206 AU, mpc 1 to 10, cpc 1 to 30.7, ly 1 to 3.26, pc 1 to 1000; speed km/h to 3,600, km/s 1 to 1,000, Mm/s 1 to 29.98, c from 0.1 c; duration µs to years as usual, then ky and My to 1000 and Gy unbounded. The km rung is 8 orders of magnitude wide (hence the scientific-km gap). The Mm/s rung runs only to 29.98 Mm/s, so a speed just under 0.1 c prints "30 Mm/s" while 0.1 c prints "0.1 c"; this follows Boss's UX.13 wording.

### Proposed ladders

Rule of thumb for every rung choice: show the largest rung whose mantissa is at least 1, as today, so mantissas are in [1, 1000) wherever neighbouring rungs are at most 1000 apart (NIST SP 811 advises a prefix that puts the number between 0.1 and 1000 [R]). Where a family has natural astronomical anchors (Earth, Jupiter, Sun) the switch points sit at the anchors' conventional ranges. All thresholds are defaults, not facts.

| Family | Base | Ladder (start threshold) | Aside / companion | Justification |
|---|---|---|---|---|
| Distance (existing, keep Boss's rungs) | m | km; AU from **1,000,000 km (0.0067 AU)**; mpc, cpc, ly, pc, kpc, Mpc, Gpc as now | parsec rungs keep "(ly)" / "(AU)" / "(km)" aside | AU is the IAU 2012 B2 exact 149,597,870,700 m [S: syrte.obspm.fr]. Mercury at 0.387 AU is "0.387 AU" in every textbook, not "5.79 × 10⁷ km". |
| Body radius | m | Stays km in scientific notation with "(x Earth radii)" (Boss, 2026-09-30) | Earth radii for planets; add Jupiter radii for gas giants | `format_body_radius_km` is a recorded decision. R_jup nominal 7.1492e7 m [S: IAU 2015 B3]. |
| Mass | kg | kg (scientific) below 0.01 lunar mass (7.3e20 kg); lunar masses to 0.1 Earth mass; Earth masses to 0.1 Jupiter mass (31.8 M⊕); Jupiter masses to 0.075 M☉ (hydrogen-burning limit, about 79 M_J); solar masses above | none needed | Prototype outputs [C]: Ceres 0.0128 lunar masses, Mercury 4.49 lunar masses, Mars 0.107 Earth masses, Neptune 17.1 Earth masses, Saturn 0.299 Jupiter masses, Proxima 0.122 solar masses, Sgr A* 4.30 × 10⁶ solar masses. Exoplanet catalogues quote both M⊕ and M_J for giants; the 0.1 M_J switch is a convention [R]. IAU 2015 B3 defines the nominal GM, not a mass, so M_sun, M_earth and M_jup depend on G. The project's M_moon (not in astropy; 7.34e22 kg [R]) should be defined as a project-nominal constant, as the IAU does. |
| Time, duration, period | s | as now, but symbols **kyr, Myr, Gyr** | none | Gy collides with the gray. Myr and Gyr are the usual astronomy and geoscience durations [S: GSA Today and stratigraphy search results; the IAU Style Manual page did not display the year rows]. The Julian year of 365.25 d is what astropy's `u.year` gives (31,557,600 s) [C]. |
| Speed | km/s | km/h; km/s from 1; Mm/s from 1,000 km/s; c from 0.1 c | none | As built (UX.13). mph only under the Customary preset (section 4). |
| Luminosity, power | W | W, kW, MW, GW, TW, PW, EW (SI prefixes); **L☉ from 1e-6 L☉** | for stars, "x L☉" stays in the aside | L_sun nominal 3.828e26 W [S: IAU 2015 B3]. Prototype [C]: Earth intrinsic 47 TW, Jupiter 335 PW, quasar 2.61 × 10¹³ L☉. Below the lowest rung use scientific text (Hawking luminosities). |
| Temperature | K | K to 3 figures; mK, µK, nK below 1 K; scientific below 1 pK | **°C and °F beside K for surface conditions only** (Boss, UX.22); stars K alone | Affine units, so the table needs an `offset`. Cross-check against `u.temperature()` (300 K gives 26.85 deg_C and 80.33 deg_F on astropy 6.0.1 and 8.0.1 [C]). |
| Pressure | Pa | nPa, µPa, mPa, Pa, kPa, MPa, GPa, TPa; scientific below nPa | **atm and psi** (shipped, PR #234); add **bar** for gas giants and **mbar** under 0.01 atm | Planetary science quotes bar, mbar, µbar (Venus about 92 bar, Mars about 6 mbar) [R], and the habitability research documents quote bar and kPa throughout (section 6). `PSI_PA = 6894.757293168361` in `constants.py` is the exact value; **do not use astropy's `imperial.psi`**, which is 6894.757388 (rounded lbf, 1.4e-8 relative error) [C]. `format_pressure_pa(2.5e-13)` currently prints "0.00000000000025 Pa" [C]. |
| Gravity (acceleration) | m/s² | m/s² primary | **g** (standard gravity 9.80665 m/s², exact) and, under Customary, ft/s² | Today only "0.38 g" is shown (`systempage._gravity_text`, `generation/planet.py`), and a bare "g" collides with the gram. Show "3.71 m/s² (0.378 g)" in detail panels, "0.378 g" in rows. Compact objects: scientific m/s². |
| Density (bulk) | kg/m³ | **g/cm³** for bodies (planet range 0.69 to 5.51 in `constants.PLANET_DENSITY`), kg/m³ for atmospheres | the other of the two | Both are already in the code (`PLANET_DENSITY` in g/cm³, `ATMOSPHERE_DENSITY` in kg/m³). |
| Number and space density | m⁻³ | particles/cm³ (nebulae, `_CLOUD_CONTENTS`), systems/ly³ (`fmt.format_density`), stars/pc³ | none | Keep the three; one number-density ladder per spatial unit. |
| Flux, heat flow | W/m² | nW/m² to kW/m² with SI prefixes; planetary heat flow in mW/m² (Earth about 87 mW/m² [R]); stellar insolation as **S⊕** (1361 W/m² nominal total solar irradiance [S: IAU 2015 B3]) | S⊕ beside W/m² | `generation/rogue.py` and `system_pages.py` show W/m²; GEN.86 and GEN.87 will add insolation. |
| Magnetic field | T | planets: nT, µT, mT (Earth 25 to 65 µT [R]); compact objects: gauss in scientific text ("1.0 × 10¹² G") | the other (1 G = 1e-4 T, checked in astropy [C]) | Pulsar and magnetar fields are quoted in gauss in the literature [R]; the DB column is `magnetic_field_gauss` and the page shows `v:.2e G`. A magnetic moment (A m²) is not stored or shown; GEN.86 stores field strength only. |
| Radiation dose | Gy (absorbed), Sv (effective) | nGy to kGy as dose per year (the habitability documents use Gy/yr, with an L_rad threshold of 10 Gy/yr); Sv/yr (mSv/yr) for human-relevant | none | astropy 6.0.1 has **no** `Gy` or `Sv` unit; 7.2.0 and 8.0.1 do [C]. The dose ladder cannot rely on astropy on Python 3.9, so define the factors by hand. Earth background about 2.4 mSv/yr and Mars surface about 0.2 to 0.3 Sv/yr (MSL/RAD) are [R]. |
| Angles | rad | degrees with no space before ° (orbital angles); arcmin, arcsec, mas for angular sizes and parallax | none | NIST SP 811: degree, minute and second signs attach to the number, the only exception to the space rule [S: NIST checklist search result]. `pitch_angle_rad` and `*_deg` columns exist. |
| Rates | kg/s, M☉/yr | "x solar masses/year" as now (quasar accretion) | none | One rung; no ladder needed. |

## 2. Libraries and the server/browser split

### astropy.units, pint, custom

| | astropy.units | pint | custom table (recommended core) |
|---|---|---|---|
| Already in the project | yes: `physics/constants.py` imports `astropy.constants` and `astropy.units`; locked | no | yes |
| Python 3.9 | astropy **6.0.1** (2024-03-26) is the last release for 3.9: 6.0.x `requires_python >=3.9`, 6.1.0 and later `>=3.10`, 7.0 `>=3.11`, 8.0 `>=3.11` [C: PyPI JSON] | 0.24.4 is the last for 3.9, but 0.24.1 fails to import in a fresh venv on 3.9 (`cannot inherit frozen dataclass from a non-frozen one`, flexparser mismatch) [C]; 0.25 needs 3.11, 0.26 needs 3.12 | any |
| numpy 2 | 6.0.1 pins `numpy<2` (so Python 3.9 stays on numpy 1.26.4, as in the lock); 6.1.x and 7.x accept numpy 2; 8.0.1 requires `numpy>=2.0` [C]. astropy 6.1.7 and 7.0.0 fail to import under numpy 2.5 (`np.in1d` removed) [C], which the lock avoids by pinning numpy 2.2.6 for 3.10. | none | none |
| Per-conversion cost, 1e5 values [C] | `Unit.to(v)` 4.5 to 5.6 µs; Quantity `.to().value` 22 to 37 µs; precomputed factor 0.1 to 0.17 µs; array `.to` 0.08 to 0.16 µs per element | 87 µs (Python 3.13 only; did not import on 3.9) | 0.17 to 0.5 µs |
| Astronomy units | au, pc, lyr, Msun, Mjup, Mearth, Rsun, Rjup, Rearth, Lsun, yr, Myr, Gyr, `imperial`; no lunar mass | some, fewer | whatever the table says |
| Temperature | `u.temperature()` equivalency (K, deg_C, imperial.deg_F) | offset units supported | `offset` field |

astropy's own formatting gives "0.999999995320789 AU", "km s⁻¹" and (with `precision=3`) "1.000 AU": no grouping, no UX.36 rule, no singular or plural. Importing `astropy.units` costs 0.6 to 0.8 s and 45 to 51 MB, already paid by the web process because `util/format.py` imports `physics/constants.py` [C].

**Recommendation:** a custom table, with astropy supplying the numbers and the cross-checks. The ladder module imports from `physics/constants.py` (which already reads astropy), and a test asserts for every rung `rung.factor == pytest.approx(u.Unit(rung.astropy_name).to(base), rel=1e-12)`. Where astropy lacks the unit (lunar mass; Gy and Sv on 6.0.1) the test skips and the factor is project-nominal. `find_equivalent_units()` is a test helper here: "every rung of the length ladder is in `u.m.find_equivalent_units(include_prefix_units=True)`" catches a typo'd symbol. `astropy.units.imperial` is fine for lb, mi and deg_F, but not for psi.

This narrows the `Library Migration Workflow.md` proposal (lines 70 to 85, 221), which lists `astropy.units` as the "unit transformation framework" that "guarantees unit correctness". What shipped in GEN.66 is the constants half (`physics/constants.py` on `astropy.constants`); the display path stays on plain floats, for the reasons in point 2 of the summary.

### Python and JS without two implementations

The page-building code is server-side (Jinja and Python), but the maps and panels (`systemmap.js`, `galaxymap3d.js`, `facilityform.js`, `phenomenonmap.js`, `galaxysystem.js`) format values live in the browser, so a JS formatter cannot be avoided. Options:

1. **Server renders everything.** Fails: hover text and live readouts (facility form slider, route lengths) change without a round trip.
2. **Intl.NumberFormat with `style: "unit"`.** Tested on Node 22 (ICU 77) [C]: supported units include meter, kilometer, mile, foot, kilogram, second to year, celsius, fahrenheit, percent, degree, kilometer-per-hour; **not** supported: kelvin, pascal, atmosphere, bar, psi, watt, au, parsec, lightyear, gauss, tesla, sievert, gray. `notation: "scientific"` gives "1.235E6", not "1.23 × 10⁶" (parts are available through `formatToParts`, so a custom renderer is possible). Compact gives "1.2M". Useful only as the digit-grouping engine, which `numberformat.js` already uses.
3. **A JS unit library** (`js-quantities` 1.8.0, `convert-units` 2.3.4, `unitmath` 1.1.1, `mathjs` 15.2.0 [C: registry.npmjs.org]). None has astronomical rungs or the UX.36 rule; the project would carry the ladder table anyway.
4. **Shared ladder table, one algorithm per language (recommended).** Python holds the table; `scripts/build_unit_ladders.py` writes `src/html/static/unitladders.js` (`export const LADDERS = {...}`, committed), and a test fails if the committed file differs from the generator's output. A generated ES module avoids `fetch` and JSON import attributes (recent browser support [R]) and fits the existing `?v=` fingerprint loading in `systemmap.js` and `components.js`. The algorithm is about 25 lines per language. Agreement is tested by a table of values run through node, which the existing `tests/test_*_format.py` already do, extended with a Hypothesis-generated batch (section 4).

Exactly this experiment was run [C]: Python and JS generic formatters consuming one JSON of the three ladders, 90,184 values (random log-uniform plus every boundary value and its neighbours). 21 mismatches, all from the half-way rounding defect (the existing `_three_figures` versus `numberformat.js` behaviour, not the prototype). After the rounding fix the twins agree.

A grid library (the TanStack Table proposal in `Web UX Development Notes.md`, section 1; see [map-ui-and-frontend-libraries.md](map-ui-and-frontend-libraries.md) for what was decided about libraries) would need no change here: its cell renderers call the JS ladder, and it sorts on the raw `data-si` value, as `datatable.js` does today with the server's `data-sort` key.

### Number rendering decisions this implies

- **Locale:** keep fixed English (comma grouping, point decimal). Boss's examples are "384,400 km" and "4.2 pc (13.7 ly)", and Markdown, wikitext and the DB must stay locale-neutral. Keep the grouping separator in one constant. NIST and the SI brochure recommend a thin space and note the comma is a decimal marker in much of the world [S, secondhand]; this is a documented deviation.
- **Rounding:** half-up on the shortest decimal form, in both languages.
- **Significant figures versus decimals:** three significant figures for ladder values (existing), but whole numbers to the unit in a rung that is not scaled (km: "384,400 km", not "384,000 km"). Keep.
- **Scientific text:** "1.23 × 10⁶" with the ISO ×, which matches ISO 80000-1 [R]. It pastes as text and spreadsheets do not parse it. The raw value goes in `data-si`; a `copy` handler that writes the plain number is possible but not worth building yet. Default: accept the limitation.

## 3. Number formatting rules (UX.36 and the rest)

**UX.36 is correctly implemented** (done in PR #457; `SCIENTIFIC_MIN_INTEGER_DIGITS = 7`, `SCIENTIFIC_MIN_DECIMAL_INTEGER_DIGITS = 5`). Verified [C]: 999,999 plain; 999,999.6 rounds to seven digits and becomes "1.00 × 10⁶"; 12,345.6 with one decimal becomes "1.23 × 10⁴"; 1,234.57 stays plain. The code docstrings, `numberformat.js` and `tests/test_number_format.py` agree with the rule. Statements that still describe the older UX.20 rule ("past 4 whole digits", "5 or more digits before the decimal point"): `docs/design/architecture.md` (the `numberformat.js` row, fixed with this note), `docs/html-interface.md` (the `numberformat.js` row), the comment above `num(x, ",.2f")` in `web/__init__.py` and the comment in `galaxystages.js` ("scientific past 4 digits (UX.20)"). The last three are code and a non-design doc, listed in the handoff.

Inconsistencies worth fixing in the same PRs:

| Issue | Evidence | Fix |
|---|---|---|
| Negative values that round to zero print "-0" in Python and JS (`format_number(-0.4)`, `-0.0`); Python `_three_figures(-0.0)` gives "-0", JS gives "0" | probe and node [C] | Normalise negative zero to "0" in both. `_whole_or_tenths` already does for temperatures. |
| Tiny values print long zero strings | summary point 4c | Scientific below mantissa 1e-3 (`_three_figures_or_tiny` already exists for atm and psi; apply it to every rung). |
| Minus sign: ASCII hyphen-minus everywhere | grep | ISO 80000-1 prefers U+2212 [R], but U+2212 breaks pasting into most spreadsheets and Python `float()`. Keep ASCII in text output; optionally render U+2212 only inside `.qty` spans, with the raw value in `data-si`. Default: keep ASCII. |
| Number and unit separated by a normal space (no NBSP anywhere in `util/format.py`, `web/lib/*.py`, `static/*.js`) | grep | A value can wrap to "384,400" / "km". NIST SP 811 requires a space and implies no break [S]. Do **not** put U+00A0 in the strings: the same strings go to Markdown and wikitext, tests compare them, and a pasted NBSP defeats `float()`. Wrap each value in `<span class="qty">` with `white-space: nowrap`. |
| Tabular numerals | `font-variant-numeric: tabular-nums` exists on `.num` and `.stat-value` in `style.css` only | Add `.qty { white-space: nowrap; font-variant-numeric: tabular-nums; }`. |
| Different units down a table column | the ladder picks per value | Add `format_column(ladder, values)` that picks one rung for the whole column (rung of the median or largest value) so rows are comparable. Sorting is by the server's `data-sort` key (`datatable.js`), unaffected. |
| Screen readers | `aria-label` is used on lists and menus only; `facilityform.js` feeds `aria-valuetext` the abbreviated "1.5 AU" | Each rung gets a `name` ("astronomical units", "light-years", "solar masses"); use it in `aria-valuetext`. For inline values, `aria-label` on a plain span is not reliably announced (ARIA gives no accessible name to generic elements) [R]; `<abbr title>` is only a tooltip. The ambiguous symbols are the ones that matter for speech: AU, ly, pc, c, g, G, K, Gyr. Default: visually hidden expansion for those only, `aria-hidden` on the abbreviation. Cost: hidden text is copied with the selection [R]. |

## 4. Unit preference and tests

**Presets (default):** `Automatic` (the ladders in section 1, the current behaviour), `Metric only` (SI base and prefixed units, no AU, ly, pc or Earth masses), `Customary` (mi, ft, lb, °F primary; mph, psi, atm). UX.42 is about wording ("Chart" instead of "Generate") and says nothing about units, and the "in-universe units" idea is not in the TODO text. A preference is a new child of UX.23.

**Where to store it:** `localStorage` key `planetgen.units`, the same pattern as `theme.js`, `bookmarks.js` and `generatefolds.js` (try/catch around every access, page works without it). Account sync later if bookmarks (UX.45) get an account store.

**How to apply it without a flash and without breaking caches:** the server always renders the default (Automatic) inside `<span class="qty" data-q="length" data-si="384400000">384,400 km</span>`; a small script reads the preference at load and re-renders each `.qty` from `data-si` with the shared JS ladder. This needs no cookie and no per-preference server render, and the page cache (`web/lib/pagecache.py`) keeps one variant. A script that runs in the head before first paint avoids layout shift for Customary users.

**Toggle patterns:** a global header menu built on `sl-dropdown` and `sl-menu-item` with checkable items (both already vendored in `components.js`; `sl-tooltip` and `sl-radio-group` are not) is the primary control. Per-value click-to-cycle costs many tab stops in tables. A tooltip showing the same value in other units ("384,400 km = 0.00257 AU = 1.28 light-seconds") serves the "contextualize" request in UX.22 without a toggle [R: pattern, not a cited source].

**Tests.** Hypothesis is already in the test extras (`hypothesis>=6.100.0`, used by `test_fuzz_*.py`; the lock caps it at 6.141.1 for Python 3.9 [C]). Properties, with the property-based tests written first:

1. Round trip: for any finite value, parsing the primary part of the output and multiplying by the rung factor is within 0.51% of the input (3 figures). 20,000 examples per ladder on length, speed and duration: pass [C].
2. Monotone rung: a larger magnitude never selects an earlier rung. 20,000 examples: pass [C].
3. Conversion identity: for every pair of rungs, `convert(convert(x, a, b), b, a)` equals x within 1e-12, and equals `(x*u.Unit(a)).to(b)` from astropy within 1e-12 where astropy has the unit.
4. Boundary table: each rung's `start`, `start*(1±1e-13)`, and 9.995, 99.95 and 999.5 times the factor (the tie-rounding cases).
5. Python/JS parity: one node call per test run on a Hypothesis-drawn list of 5,000 values per ladder; compare strings. The existing tests already shell out to node and skip when it is missing.
6. Never raises and never prints "nan" or "inf" (the current `dash_unless_finite` behaviour, TEST.53).
7. Generated-file drift: `unitladders.js` equals the generator's output.

## 5. Recommended module API and migration order

### Data structure

```python
# src/planetgen/util/unitladder.py  (Python 3.9 compatible)
from dataclasses import dataclass
from typing import Optional, Tuple

@dataclass(frozen=True)
class Rung:
    symbol: str                      # "AU"
    factor: float                    # base units per one of this unit (from constants.py / astropy)
    start: Optional[float] = None    # smallest |value| in base units shown in this rung; default = factor
    singular: Optional[str] = None   # label when the shown number is exactly 1 ("year")
    name: Optional[str] = None       # spelled out, for aria and tooltips ("astronomical units")
    offset: float = 0.0              # temperature only: K -> deg C
    aside: Tuple[Tuple[float, str], ...] = ()   # ((min |value| in base, rung symbol), ...) high to low

@dataclass(frozen=True)
class Ladder:
    quantity: str                    # "length"
    base: str                        # "m"
    rungs: Tuple[Rung, ...]          # smallest first
    companions: Tuple[str, ...] = () # always-shown secondaries, e.g. ("atm", "psi") for pressure
```

### Functions

```python
def format_quantity(ladder, value, *, system="auto", companions=True, empty="–") -> str
def quantity_parts(ladder, value, *, system="auto") -> Parts   # Parts(number, symbol, name, raw_si, text)
def quantity_html(ladder, value, *, system="auto") -> str      # <span class="qty" data-q=... data-si=...>
def format_column(ladder, values, *, system="auto") -> list[str]   # one rung for the column
def pick_rung(ladder, size, system="auto") -> Rung
def convert(value, from_symbol, to_symbol, ladder) -> float
```

`format.py` keeps its public names as one-line wrappers (`format_distance_m = partial(format_quantity, LENGTH)`) so the many call sites (about 40 files use them) do not change in step 0. JS: `unitladder.js` exports the same functions over `LADDERS` from the generated `unitladders.js` and re-exports the old names (`formatDistanceKm`, `formatSpeedKms`, `formatPeriodYears`, `LIGHTYEAR_M`, `PARSEC_M`) so `galaxymap3d.js`, `facilityform.js`, `systemmap.js`, `phenomenonmap.js`, `galaxysystem.js` and `generatebuttons.js` keep working.

Example outputs (prototype) [C]:

```
format_quantity(LENGTH, 6.4e12)      -> "42.8 AU"
format_quantity(LENGTH, 3.9e16)      -> "1.26 pc (4.12 ly)"
format_quantity(MASS, 5.683e26)      -> "0.299 Jupiter masses"
format_quantity(PERIOD, 3.1557600e7) -> "1 year"
```

### Migration order

One quantity family per PR; every step keeps pages identical except where stated. Function and file names are given instead of line numbers.

0. **Module and generated JS, no output change.** New: `src/planetgen/util/unitladder.py`, `src/planetgen/util/ladders.py` (data), `scripts/build_unit_ladders.py`, `src/html/static/unitladder.js`, generated `src/html/static/unitladders.js`, `src/tests/test_unit_ladders.py`. Changed: `util/format.py` (wrappers), `html/static/distance.js`, `speed.js`, `period.js` (re-export), `physics/constants.py` (pin the constant sets; the first import of astropy must follow the `set()` calls, so also check `planetgen/__init__.py` and any module that imports astropy earlier), `tests/test_astropy_constants.py` (tighten `HYDROGEN_ATOM_MASS_KG` to 1e-12 and add a "constant set is CODATA 2018 / IAU 2015" assertion). Include the rounding fix (summary point 4a) here, with the test table, since it is a deliberate output change in a handful of tie cases.
1. **Distance:** AU start at 1,000,000 km; `systemview3d.js` and `galaxysystem.js` orbit fields onto the ladder (radii stay scientific km); the local wrapper in `systemmap.js` removed; `searchpage.py`, `cli/query.py` and the `population/facilities.py` messages (they print raw km).
2. **Time:** kyr/Myr/Gyr symbols; `system_pages.py` age rows (`age_years`), `generation/star.py` age text, `format.format_age_string`, `population_pages.py`, and the period fields in `galaxysystem.js` and `systemview3d.js`.
3. **Temperature, pressure, gravity (UX.22 remainder):** replace `system_pages._temperature_text` and `_pressure_text` (rogue planets: K and °C only, bar only; the shipped `format_temperature_k` and `format_pressure_pa` exist but this page does not use them), `systempage._gravity_text`, the gravity line in `generation/planet.py`, `sector_page.py` (`int(temperature_k) K`), `star.py` (temperature text). Pressure gains sub-Pa rungs and bar/mbar.
4. **Mass:** `tabledisplay.format_star_mass`, `format_body_mass`, `format.format_relative_to_sol`, and the black-hole, neutron-star, rogue and quasar mass rows in `system_pages.py`.
5. **Luminosity, power, flux:** `tabledisplay.format_star_luminosity`, the luminosity rows in `system_pages.py`, `internal_heat_flux_w_m2`, `generation/rogue.py`, the L☉ text in `galaxymap3d.js`.
6. **Density, number density, magnetic field, angles:** `fmt.format_density`, `_CLOUD_CONTENTS`, `magnetic_field_gauss`, the density and field text in `galaxymap3d.js`. GEN.86 and GEN.87 add the dose and field ladders at birth.
7. **JS stragglers:** the 17 `toFixed`/`toPrecision` calls in `galaxymap3d.js`, `galaxystages.js`, `galaxystageview.js`; a lint-style test: no `toPrecision` or `toLocaleString` on a measurement outside `numberformat.js`.
8. **Preference menu and `.qty` re-render** (section 4); then UX.30's property grids call `quantity_html`.

UX.30 (structured HTML in place of the Markdown render) should wait for step 0, so the new property grids emit `.qty` spans from the start and do not re-hard-code strings. UX.32 touches only the row chips (distance, period, gravity), which come from steps 0 and 3.

## 6. Units the research documents use

Boss's research documents in this folder carry no formatting rules, but they fix which units the habitability model works in, and the display ladders should agree with them. `Web UX Development Notes.md` has no unit or number-formatting guidance (it proposes libraries for grids, 3D tiles, logs, progress bars and queues); nothing in it needed carrying over beyond the grid-library point in section 2.

| Source document | Units it uses | Consequence for the ladders |
|---|---|---|
| `Atmospheric Toxicity.md`, `Chemical Habitability.md` | partial and total pressures in bar (0.005 to 100+ bar), with kPa beside it for the human range (50 to 250 kPa, `Atmospheric Toxicity.md`) | The pressure ladder needs bar and mbar rungs (section 1), and partial pressures (pO₂, pCO₂) should use the same ladder as total pressure. |
| `Naturally Occuring Ionizing Radiation.md` | atmospheric column mass in g/cm² (X, with 1000 g/cm² at about 1 bar), human-relevant doses in mSv/yr, field in µT (a 0.05 bar, 30 µT case) | Column mass is a new quantity (g/cm², or kg/m²); doses split between Sv/yr for people and Gy/yr for the habitability thresholds. |
| `Planetary Habitability Index.md` | unshielded dose above 10 Gy/yr as the harshest class, 50 mSv/yr as the occupational limit | The dose ladder shows Gy/yr and Sv/yr as separate quantities, never a bare "Gy" for a duration (point 4d). |
| `Planetary Habitability and Speculative Xenobiology.md` | radiolytic energy deposition in eV g⁻¹ s⁻¹ | Internal to the habitability score; it needs no ladder unless a page shows it. |

[habitability-index.md](habitability-index.md) already flags that the L_rad threshold (10 Gy/yr) disagrees with the Mars figure in the source table; that is a unit-and-magnitude question for the dose ladder when GEN.87 lands, and is not settled here.

## Sources

- IAU 2015 Resolution B3 (nominal solar and planetary conversion constants): https://www.iau.org/common/Uploaded%20files/IAUGA2015-Resolution-B3-recommended-nominal-conversion.pdf , https://www.pas.rochester.edu/~emamajek/IAU/IAUres_B3.pdf , https://arxiv.org/pdf/1605.09788 (search results only: R_sun 6.957e8 m, GM_sun 1.3271244e20, L_sun 3.828e26 W, S 1361 W/m², R_jup(eq) 7.1492e7 m, "nominal values are exact conversion factors, not best estimates")
- IAU 2012 Resolution B2 (au = 149,597,870,700 m, symbol "au"): https://syrte.obspm.fr/IAU_resolutions/IAUResol_2012_0.html , https://syrte.obspm.fr/~capitaine/AU.html
- NIST SP 811 chapter 7 and checklist (space between number and unit, plane-angle exception, comma not for grouping): https://nist.gov/pml/special-publication-811/nist-guide-si-chapter-7-rules-and-style-conventions-expressing-values , https://www.nist.gov/document/sp-811-2008-checklist-reviewing-manuscriptspdf (snippets via search)
- SI brochure digit-grouping excerpt (8th edition, via mailing-list quote in search results): https://febo.com/pipermail/time-nuts/2009-October/041350.html
- IAU Style Manual (Wilkins 1989) reprint: https://iauarchive.eso.org/publications/proceedings_rules/units/ ; Myr/Gyr usage: https://rock.geosociety.org/gsatoday/archive/22/2/pdf/i1052-5173-22-2-28.pdf , https://www.micropress.org/microaccess/stratigraphy/issue-260/article-1641
- Intl.NumberFormat: https://github.com/tc39/proposal-unified-intl-numberformat/blob/master/README.md , https://v8.dev/features/intl-numberformat (plus Node 22 / ICU 77 tests)
- PyPI JSON: https://pypi.org/pypi/astropy/json , https://pypi.org/pypi/pint/json , https://pypi.org/pypi/hypothesis/json ; npm registry for js-quantities, convert-units, unitmath, mathjs

## Evidence notes

[C] results came from throwaway scripts in the research sandbox, not in the repository: equivalence of the data-driven ladders with the existing formatters (600,000 values, 0 mismatches); the Python/JS cross-check (90,184 values, 21 mismatches, all tie rounding); tie-rounding tests (`Decimal(repr(x))` half-up equals Intl on 400,000 values); astropy constants per version; missing `Gy` and `Sv` in astropy 6.0.1; the `imperial.psi` value; benchmarks (Python 3.9 + astropy 6.0.1 + numpy 1.26.4, and Python 3.13 + astropy 8.0.1 + numpy 2.5.3); Intl unit support on Node 22; the Hypothesis properties; the pint 0.24.1 import failure on 3.9; astropy 6.1.7 and 7.0.0 failing under numpy 2.5. The repository checks (UX.36 limits, `PSI_PA`, the lock's astropy versions, the 189-line JS twins, the `HYDROGEN_ATOM_MASS_KG` tolerance, the stale "4 whole digits" statements) were re-verified against the checkout on 2026-10-09.

[R], to confirm when web access allows: NIST SP 811 "prefix so the number lies between 0.1 and 1000"; ISO 80000-1 minus sign, thin-space grouping and the "×" sign; IAU 2015 Resolution B2 on the parsec (astropy's reference string for `pc` reads "Derived from au + IAU 2015 Resolution B ...", truncated at 40 characters in the output, which supports it); Moon mass 7.34e22 kg; Earth heat flow 87 mW/m²; Earth field 25 to 65 µT; Earth and Mars dose figures; exoplanet-catalogue use of M⊕ versus M_J and the 0.1 M_J switch; hydrogen-burning limit 0.075 M☉; planetary pressures in bar and mbar; ARIA naming of generic elements and the copy behaviour of visually hidden text; browser support for JSON import attributes; unit-toggle UX patterns; the Python 3.9 end-of-life date; whether astropy 6.1.x and 7.0 have `Gy` and `Sv` (not tested because they fail to import under numpy 2.5).
