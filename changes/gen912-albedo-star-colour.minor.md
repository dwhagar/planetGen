**Headline:** The Bond albedo of an ice-covered or ocean world now follows the colour of its star, so planets of K and M dwarfs are a little warmer at the same flux.

### Changed
- **The colour of the star shifts the albedo of a bright world (GEN.912).** Ice, snow and ocean under cloud reflect the blue light of a hot star well and the near-infrared light of a cool one badly, so the Bond albedo range of a frozen world is now scaled by 0.79 (M dwarf, 3,410 K), 0.935 (K), 1.0 (Sun-like) and 1.05 (F, 6,900 K), and a wet world takes half of the ocean-row shift of Shields et al. 2013 (Table 2); the range of a dry world is unchanged. The factors are registered stand-ins (SAB1 to SAB3) from Shields et al. 2013 and Joshi and Haberle 2012. A frozen planet of an M dwarf comes out about 20 percent more absorbing than before and so a few kelvin warmer at the same flux; a Sun-like one does not change.
- Generator epoch bump (numbered at merge, `--affects fill`): the albedo, temperature, climate and everything derived from them change for bodies of K and M dwarfs (and F stars).

generation-output: changed
