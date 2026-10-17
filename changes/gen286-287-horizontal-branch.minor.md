**Headline:** Old stars now have a horizontal branch and an asymptotic giant branch, and their red giant branch spends its time at each luminosity as a real one does, so the bright-giant counts match real ones.

### Added
- **A horizontal branch and an asymptotic giant branch (GEN.286).** A star lighter than 2.2 solar masses burns helium in its core for about a hundredth of its main-sequence life as a class III star of 40 to 60 solar luminosities, then climbs the asymptotic giant branch for a five-hundredth of it, before it becomes a white dwarf. The horizontal-branch colour (red clump, RR Lyrae, blue or extreme, by [Fe/H] and age) sets its temperature, and the star keeps the [Fe/H] its colour was drawn from. They show in the scatter at a floor of 50 to 60 solar luminosities and in every sector fill (halo, bulge, clusters), with "Horizontal Branch Giant" or "Asymptotic Giant" in the type text. The pre-placement tables and the mass bands mirror the new phases.

### Changed
- **The red giant branch runs from its base to its tip with the real luminosity function (GEN.287).** The time a giant spends per unit of ln L falls as L^-0.9, from its base (about 3 solar luminosities for an old turn-off star) to 2,500, instead of an even spread over 50 to 1,000 with a 2% tip. Most giants are now dim K giants and the bright ones are far fewer: the halo scatter at a floor of 500 solar luminosities stores about 190,000 rows, not 3.8 million, and its counts agree with the evolutionary flux of an old population. A galaxy scattered before this has too many bright giants: `planetgen plan --redo-scatters mass luminosity` draws them again.
- **A star that ends its giant phases becomes a white dwarf later, by 1.2% of its main-sequence life** (the horizontal branch and the asymptotic giant branch lie between), and a subgiant brightens to 3.5 times its main-sequence luminosity, not twice.

TODO items resolved: GEN.286, GEN.287. Added: none.
