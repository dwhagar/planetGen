# The isotope record of a body (GEN.1377, GEN.1420)

All real materials exist (Boss, 2026-11-08): deuterium is a share of hydrogen, not a fictional substance, and every stable or long-lived isotope of an element a body holds can be drawn. This note says what a body stores for it, which laws move the shares away from the natural mix, and what is left for the second part. The code is `physics/isotope_record.py`; the stand-ins are section 6.233 of the gap-filling plan (ISO1 to ISO9).

## 1. What is stored

A planet, moon or rogue planet carries a list of rows, `isotopes`, kept in `planet_isotopes`, `moon_isotopes` and `rogue_planet_isotopes` (schema v235): `element`, `isotope`, `fraction` (the atom fraction of the element; the rows of one element add up to 1) and `origin` (`natural`, `physical`, `decay`, `radiogenic`). It is built after the air, the interior and the ground, on the body's own streams (no draw of the run's moves), and again when an edit moves the air or the interior. A body saved before v235 has no rows.

## 2. Which elements a body holds

The elements of what it is made of, read from the stored composition: the rock (O, Mg, Si, Al, Ca, Fe and each minor oxide with a share above zero; C and Si for a carbon rock; K, Th and U where the body stores a radiogenic content), the core (Fe and each light element with a share), the water (H, O), each ice (the elements of its formula), an envelope (H and He, and C, N, O where it holds heavy elements) and every gas of the air.

## 3. The laws

| Element | Law | Sources |
|---|---|---|
| H | D/H of the star's birth gas = primordial 2.527e-5 times `exp(-A 10^[Fe/H])`, A fixed by the protosolar 2.1e-5; ice and water enriched by `exp(T0 (1/T - 1/182 K))`, fitted to the ocean (1.56e-4 at 182 K) and comet 67P (5.3e-4 at 30 K); the body's D/H is the hydrogen-atom-weighted mean of its envelope, water and ices, with one log-normal scatter | [R0088, R5171, R5173, R5174] |
| C | 12C/13C from the gradient 6.21 D + 18.71 read at the radius whose gas has the star's metallicity, moved between the local interstellar 68 and the Solar System's 89 by the star's age | [R5175] |
| N | 14N/15N from the gradient 21.1 D + 123.8, between the local ~400 and the protosolar 441, then the solids raised in 15N by the same kind of law as D/H (fitted to the terrestrial 274 and the comets' ~140); an envelope keeps the gas value | [R5172, R5177] |
| O | the table's natural mix (no source read for the oxygen gradient; the mass-independent anomalies are left out) | |
| Ar | 36Ar/38Ar = 5.501 (primordial); 40Ar = the 40K that decayed times 10.72 percent times the interior's degassed share, over the moles of argon in the air, at most 99.96 percent | [R5176] |
| K, U | the natural mix decayed back to the Solar System's birth (a stand-in for every star) and forward to the body's age | `radiogenic.formation_mix` |
| all others | the table's natural mix (NIST and IUPAC) | |

## 4. What a generated sample shows

A solar-metal planet 4.5 Gyr old at 1 AU draws D/H near 2e-4 (the ocean's 1.56e-4 scattered), 13C/12C 1/89, 15N/14N 1/274, a 99.96 percent 40Ar air and 235U 0.76 percent of its uranium. An icy body at 30 K draws D/H near 5e-4. A hydrogen giant keeps the birth gas's 2.1e-5. A metal-rich star's bodies are richer in 13C and 15N, a metal-poor star's in deuterium.

## 5. Limits, and the second part

The short-lived isotopes (Fe-60, Mn-53, I-129, Pu-244), the xenon excess they leave, the cosmogenic tritium, C-14, Be-10 and Al-26 (GEN.1421), the daughters of the long-lived chains (lead, strontium, neodymium), the mass-independent oxygen anomalies, the split between a body's reservoirs (an ocean against a mantle), and belts, comets and stars are not in this part. The page lists the ratios the record moves (D/H, 13C/12C, 15N/14N, 40Ar/36Ar, 235U/238U, 40K/39K) and the count of isotopes recorded.
