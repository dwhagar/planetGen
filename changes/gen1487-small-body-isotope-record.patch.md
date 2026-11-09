**Headline:** Belts, comets, interstellar comets and interstellar fields now show a full isotope record on their property sheet and in the API.

### Changed

- A belt, a comet, an interstellar comet and an interstellar field list D/H, 13C/12C and 15N/14N (where the ice holds the element) and the count of isotopes, worked out from their stored ice share, the host star's [Fe/H] and age and the temperature the ice formed at, through the same record and ice fractionation law as a planet's (GEN.1487). Before, a belt or comet showed only two ice lines.
- Nothing is stored and nothing generated changes: the record is derived each time from stored values and a seeded stream named by the body, so a body always shows the same line.
- New stand-ins ISO18 to ISO20 (the ices a small body holds, the budget margin at the water line, the formation temperature and metallicity of a body with no star) are in section 6.236 of the registry.

generation-output: unchanged
generation-affects: none
