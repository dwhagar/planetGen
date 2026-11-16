**Headline:** Ordinary red dwarfs are no longer labelled Subdwarf; only metal-poor stars ([Fe/H] -0.5 or lower) can be.

### Fixed

- GEN.1635: `stellar_evolution.yerkes_class_of` called a main-sequence star Subdwarf (class VI) when it lay 0.3 dex under the dwarf sequence's luminosity at its temperature. The model makes M dwarfs under-luminous from [Fe/H] -0.1 on (0.3 to 0.5 dex at -0.2), so 13 of 270 ordinary single M3 to M5 stars of [Fe/H] -0.08 to -0.34 were labelled Subdwarf, while real subdwarfs are metal-poor dwarfs (sdM from [Fe/H] about -0.5, extreme ones near -1 and lower). A star is now VI only when its birth [Fe/H] is also at or below the new `tuning.SUBDWARF_MAX_FE_H` (-0.5). In a draw of 3,000 M dwarfs of 0.1 to 0.45 solar masses at the Sun's place, 301 were labelled Subdwarf before and none now. Forced subdwarf types (halo [Fe/H] -2 to -1) are unchanged. Stand-in PS14 now says so.

generation-output: changed
generation-affects: scatter, fill
