**Headline:** The smallest galaxy the Scale dial offers, 281 solar masses of stars, is now as small as the real one it is half of (Segue 1), not nine times too large.

### Changed
- GEN.612: the thin disk's scale length came from [R4048]'s late-type size-mass fit down to every mass, but that fit starts at about 1e8.8 solar masses of stars: at the Scale dial's floor of 281 solar masses it gave a half-light radius of 220 pc, nine times Segue 1's 23.9 pc. Under 10^8.5 solar masses (`shape_draw.DWARF_SIZE_EDGE_SOLAR`) the median is now the field dwarfs' law of [R0655] (`log r_e = 0.667 + 0.296 log M*`, fitted over log M* 5.5 to 8.5, scaled by 1.018 so it meets [R4048]'s at the edge): 25 pc at 281 solar masses, 558 pc at 1e7, 1.55 kpc at the edge. Below log M* 5.5 the law is continued as a labelled placeholder (rows DG1 to DG3 of the new 6.332; Segue 1 is 0.11 dex under it) and a dwarf's height is the thin disk's flatness (7.3) until the dwarf families are built. The halo needs no extension: a galaxy under 4.6e6 solar masses draws none (Boss, 7:33 PM EDT on 2026-10-21) and the relation of the mass chain stands above it (row HF3). A test sweeps the dial's whole range, 281 to 1.55e11, for finite, ordered laws. Galaxies of 3.2e8 solar masses and more, the Auto window included, draw as before.

generation-output: changed
generation-affects: plan
