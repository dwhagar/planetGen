**Headline:** Every planet, moon and rogue planet now carries an isotope record: the shares of each isotope of each element it holds, with deuterium, carbon-13, nitrogen-15, argon-40 and the decaying potassium and uranium following physical laws, and the record shown on the property sheet (GEN.1377, GEN.1420 part 1).

### Changed

- New `physics/isotope_record.py` builds the record from the body's stored composition on its own streams: the table's natural mix by default and, for hydrogen, carbon, nitrogen and argon, a law of the star's metals and age and of the temperature the ice formed at (D/H from the primordial value less the deuterium star formation burned, enriched in ice as the temperature falls, fitted to the ocean and comet 67P; 12C/13C and 14N/15N from the galactic gradients; 36Ar/38Ar from the solar wind; 40Ar from the potassium that decayed and the interior let out). Potassium and uranium isotopes are those of the body's age.
- Schema v235 (migration 0235): the child tables `planet_isotopes`, `moon_isotopes` and `rogue_planet_isotopes`. A body saved before it has no record and nothing is backfilled.
- The property sheet of a planet or moon lists the ratios the record moves and the count of isotopes. `bulk_draw.formation_temperature_k` is now public.
- Stand-in section 6.236 (ISO1 to ISO9) in the registry, design note `composition/isotope-record.md`, citations R5173 to R5179.
- The short-lived isotopes and the cosmogenic ones (GEN.1377, GEN.1421) are the second part.

generation-output: changed
