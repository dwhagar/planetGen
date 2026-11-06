**Headline:** A build no longer dies on a rock whose dayside is hotter than 6000 K.

### Changed

- The vapor over a molten surface (`physics/mixture_lava.py`) is worked out at no more than `MELT_MAX_K`, 6000 K, the top of NASA Glenn's fits of the liquid oxides and of 47 of the 84 gases (GEN.1294). A build had died once with `SiO2(L): no thermodynamic fit at 6881 K (the fits cover 1996 to 6000 K)` for a hot-zone rock whose dayside (`tp_k`) or a white-dwarf planet's exosphere (`ts_k`) was over 6000 K; a rerun passed because it drew a cooler body. The mixture still reports the surface's own temperature, and `species_thermo` still raises for any other call outside its fits. The clamp is a stand-in, row AM95 of the registry (tier 6, tweak level 0).
- Every surface at or under 6000 K gives the same vapor as before, so no stored value changes; only a body that used to stop the build now has a vapor.

generation-output: unchanged
