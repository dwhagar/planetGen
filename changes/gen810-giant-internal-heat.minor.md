**Headline:** A giant that has rained helium into its depths now radiates more heat than one that has not, and the planet stores it.

### Changed
- **A giant's internal heat reads its interior (GEN.810).** A gas giant, ice giant or brown dwarf stored no heat of its own. It now stores `internal_heat_flux_w_m2`, the heat it radiates per square metre: Jupiter's 100 K of own heat scaled by the square root of the mass, plus, for a giant whose helium rained out of the molecular envelope, the energy the helium gave up sinking (read back from the stored layers, `helium_depletion` and `has_helium_shell`; half of it leaves as heat, spread over the star's age). A Saturn-like giant with a helium shell radiates about 1.7 W/m2 against 0.5 without the rain (1.88 measured); a Jupiter-like giant a few percent more than its base. The one-bar temperature that decided the rain is unchanged.
- **Stand-ins recorded.** Rows GH1 to GH4 in the stand-in registry (GEN.725): the base law, the share of the energy that leaves as heat, the age it is spread over and the approximations.
generation-output: changed
