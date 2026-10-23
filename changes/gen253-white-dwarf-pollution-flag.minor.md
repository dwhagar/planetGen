**Headline:** A cool white dwarf now carries a pollution flag at the observed share, drawn from the star's own stream so nothing is stored.

### Added

- `white_dwarf_polluted` and `StarMakeup.detail["polluted"]` (GEN.253): a white dwarf cooler than 10,000 K is polluted by accreted rocky debris with probability 0.375, the middle of the 25 to 50 percent observed; a hotter one never is. A new last draw, `StarDraws.u_pollution`, decides it, so no earlier draw moves.
- `multistar.draw_white_dwarf_pollution(teff_k, key)`, the generator the `white-dwarf-pollution` row of the system-type registry names; it reads the registry's default share and agrees with the makeup.
- The "cool white dwarf presets" of the item: a cool white dwarf is already a forced type (`--star-type A5VII` and any class from A to M with class VII), and the Generate page and command line list no star types to add to, so no preset was added.
