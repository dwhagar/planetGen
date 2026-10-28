**Headline:** The magnetopause distance on a world's page is now in radii of that world, as its label says, so Jupiter and every world that is not Earth-sized read right.

### Changed
- **Magnetopause in planet radii (GEN.843).** `magnetopause_rp` was stored in Earth radii but labelled planet radii (the unit registry's `R_planet` and the sheet's "Magnetopause distance"), so any body not Earth-sized read wrong by the ratio of radii. `magnetism.magnetopause_rp` now divides the standoff by the body's own radius, which agrees with the metres the radiation belts work from. A test reads the stored value for Earth, Mars and Jupiter back against those metres.

TODO items resolved: GEN.843. Added: none.
