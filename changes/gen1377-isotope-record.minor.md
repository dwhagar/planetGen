**Headline:** Every planet, moon and rogue planet now carries an isotope record: the shares of each isotope of each element it holds, with deuterium, carbon-13, nitrogen-15, argon-40 and the decaying potassium and uranium following physical laws, and the record shown on the property sheet (GEN.1377, GEN.1420, GEN.1421): also the short-lived isotopes of the star's birth gas and the tritium and radiocarbon the cosmic rays make in a body's air.

### Changed

- New `physics/isotope_record.py` builds the record from the body's stored composition on its own streams: the table's natural mix by default and, for hydrogen, carbon, nitrogen and argon, a law of the star's metals and age and of the temperature the ice formed at (D/H from the primordial value less the deuterium star formation burned, enriched in ice as the temperature falls, fitted to the ocean and comet 67P; 12C/13C and 14N/15N from the galactic gradients; 36Ar/38Ar from the solar wind; 40Ar from the potassium that decayed and the interior let out). Potassium and uranium isotopes are those of the body's age.
- Schema v235 (migration 0235): the child tables `planet_isotopes`, `moon_isotopes` and `rogue_planet_isotopes`. A body saved before it has no record and nothing is backfilled.
- The property sheet of a planet or moon lists the ratios the record moves and the count of isotopes. `bulk_draw.formation_temperature_k` is now public.
- Stand-in section 6.236 (ISO1 to ISO17) in the registry, design note `composition/isotope-record.md`, citations R5173 to R5182.
- A belt's and a comet's property sheet list the D/H and 15N/14N of their ices, worked out from the stored ice share, the star's [Fe/H] and age and the temperature the ice formed at (a belt's from its distance, a comet's the cold anchor), nothing stored.
- Rows of origin `extinct` (26Al, 60Fe, 53Mn, 129I, 244Pu: the star's birth ratio and what is left of it) and `cosmogenic` (3H and 14C: the steady-state share in the air and ocean from the production of cosmic rays, the air column, the targets and the field). 10Be and 26Al made in the air wait for a source.

generation-output: changed
