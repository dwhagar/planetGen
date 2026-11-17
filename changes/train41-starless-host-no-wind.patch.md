**Headline:** A body with no star has no wind again, so a starless host's planet has no magnetopause standoff, and the subdwarf [Fe/H] cut is a signed constant.

### Fixed

- GEN.1625 follow-up: `magnetism.wind_pressure_pa` held every wind at the interstellar medium's pressure, so a rogue planet's starless host (luminosity 0) got a wind and `magnetopause_rp` a value where it was `None` (GEN.1111). A star with no light (a starless host, a black hole nothing feeds) now gives no wind; the floor stays for a quiet star and for a white dwarf. `test_fuzz_utils_and_config` lists `SUBDWARF_MAX_FE_H` (GEN.1635, -0.5) among the constants that are negative by nature.

generation-output: unchanged
