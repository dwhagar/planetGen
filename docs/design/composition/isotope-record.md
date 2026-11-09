# The isotope record of a body (GEN.1377, GEN.1420)

All real materials exist (Boss, 2026-11-08): deuterium is a share of hydrogen, not a fictional substance, and every stable or long-lived isotope of an element a body holds can be drawn. This note says what a body stores for it, which laws move the shares away from the natural mix, and what is left for the second part. The code is `physics/isotope_record.py`; the stand-ins are section 6.236 of the gap-filling plan (ISO1 to ISO9).

## 1. What is stored

A planet, moon or rogue planet carries a list of rows, `isotopes`, kept in `planet_isotopes`, `moon_isotopes` and `rogue_planet_isotopes` (schema v235): `element`, `isotope`, `fraction` (the atom fraction of the element; the rows of one element add up to 1) and `origin` (`natural`, `physical`, `decay`, `radiogenic`). It is built after the air, the interior and the ground, on the body's own streams (no draw of the run's moves), and again when an edit moves the air or the interior. A body saved before v235 has no rows.

## 2. Which elements a body holds

The elements of what it is made of, read from the stored composition: the rock (O, Mg, Si, Al, Ca, Fe and each minor oxide with a share above zero; C and Si for a carbon rock; K, Th and U where the body stores a radiogenic content), the core (Fe and each light element with a share), the water (H, O), each ice (the elements of its formula), an envelope (H and He, and C, N, O where it holds heavy elements) and every gas of the air.

## 3. The laws

| Element | Law | Sources |
|---|---|---|
| H | D/H of the star's birth gas = primordial 2.527e-5 times `exp(-A 10^[Fe/H])`, A fixed by the protosolar 2.1e-5; ice and water enriched by `exp(T0 (1/T - 1/182 K))`, fitted to the ocean (1.56e-4 at 182 K) and comet 67P (5.3e-4 at 30 K); the body's D/H is the hydrogen-atom-weighted mean of its envelope, water and ices, with one log-normal scatter | [R0088, R5173, R5175, R5176] |
| C | 12C/13C from the gradient 6.21 D + 18.71 read at the radius whose gas has the star's metallicity, moved between the local interstellar 68 and the Solar System's 89 by the star's age | [R5177] |
| N | 14N/15N from the gradient 21.1 D + 123.8, between the local ~400 and the protosolar 441, then the solids raised in 15N by the same kind of law as D/H (fitted to the terrestrial 274 and the comets' ~140); an envelope keeps the gas value | [R5174, R5179] |
| O | the table's natural mix (no source read for the oxygen gradient; the mass-independent anomalies are left out) | |
| Ar | 36Ar/38Ar = 5.501 (primordial); 40Ar = the 40K that decayed times 10.72 percent times the interior's degassed share, over the moles of argon in the air, at most 99.96 percent | [R5178] |
| K, U | the natural mix decayed back to the Solar System's birth (a stand-in for every star) and forward to the body's age | `radiogenic.formation_mix` |
| all others | the table's natural mix (NIST and IUPAC) | |

## 4. What a generated sample shows

A solar-metal planet 4.5 Gyr old at 1 AU draws D/H near 2e-4 (the ocean's 1.56e-4 scattered), 13C/12C 1/89, 15N/14N 1/274, a 99.96 percent 40Ar air and 235U 0.76 percent of its uranium. An icy body at 30 K draws D/H near 5e-4. A hydrogen giant keeps the birth gas's 2.1e-5. A metal-rich star's bodies are richer in 13C and 15N, a metal-poor star's in deuterium.

## 5. The short-lived and cosmogenic isotopes (GEN.1377, GEN.1421)

**Extinct tracers.** A body with rock carries 26Al, 60Fe, 53Mn, 129I and 244Pu as rows of origin `extinct`: `initial_ratio` is the ratio to 27Al, 56Fe, 55Mn, 127I and (for plutonium, which has no stable isotope) 238U in the gas the star was born from, and `fraction` is what is left per atom of that reference element at the star's age, a number that is nothing for any star older than a few hundred million years. 26Al is the ratio `radiogenic` already draws. 60Fe is a sum of decayed additions about the gamma-ray mean of the galaxy, built the way 26Al is (`birth_radioisotopes`). 53Mn, 129I and 244Pu use the same sampling about a tenth of the early Solar System's ratios (`ISO10`, `ISO11`; mean lives from the review's Table 2). Every body of a star shares the star's ratios.

**Cosmogenic tracers.** A body with air carries 3H and 14C as rows of origin `cosmogenic`: `fraction` is the steady-state atoms per atom of the element in the pool the isotope mixes into (the air, and for hydrogen the surface ocean's water, for carbon the air and ocean together at the Earth's ratio). The production is the Earth's global rate times the share of the cascade the air column catches (`1 - exp(-X/150 g/cm2)`), the targets per gram of air (tritium from N and O, radiocarbon from N alone) and the field's shielding (a body with no dipole gets 2.09 times the Earth's rate). The Earth's own air gives 1.6 and 0.345 atoms per cm2 per second and 14C/C of 1.3e-12, 3H/H of 1e-20. The sheet lists 3H/1H and 14C/12C beside the others.

**Small bodies** (derived, not stored; GEN.1487): a belt, a comet, an interstellar comet and an interstellar field have no record table. `isotope_record.small_body_record` works one out from the values the row stores (its ice share, the host star's [Fe/H] and age, the temperature the ice formed at) through `record_of`, the planet laws of section 3 with no second law, so the sheet and the API list the same ratio lines (D/H, 13C/12C, 15N/14N, the count) and a census can call the function on any row it draws. The ice is split into water, NH3, CO2, CH4, N2 and CO by the disk budget at the formation temperature (`ISO18`, `ISO19`), so a warm belt holds water only and shows D/H alone, and nitrogen lines appear where nitrogen ice froze. A belt's temperature comes from its distance and its snow line, a comet's and an interstellar comet's is the 30 K of the cold anchor (`ISO17`, `ISO20`), a field's the water line's; a body with no host star takes the gas [Fe/H] of its galactic radius and the Earth's reference age (`ISO20`). The scatter is one draw named by the body (system and distance, or its id), so a body always shows the same line. A member of a belt or field is given its own rock, metal, water and ice fractions.

## 6. Limits, and what is next

10Be and 26Al made in the air (no global rate in hand; beryllium has no pool the record holds), cosmic-ray production in an airless or ocean-covered surface, the stellar wind and the star's own flares (the Earth's galactic flux is used for every star), the xenon excess that 129I and 244Pu leave, the daughters of the long-lived chains (lead, strontium, neodymium), the mass-independent oxygen anomalies, the split between a body's reservoirs (an ocean against a mantle), star surfaces and the sheets of the individual members of a belt or field are not in this record yet. The page lists the ratios the record moves away from the natural mix.
