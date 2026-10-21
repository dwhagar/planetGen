**Headline:** A galaxy's neutron stars are now what its massive stars leave when they do not become black holes, instead of a fixed number per star.

### Changed
- The neutron star density at a place is the black hole sum with the survivors of the failed-supernova share (`1 - F`) and the galaxy's bound share, so a galaxy whose massive stars fail more makes more black holes and fewer neutron stars, and a metal-rich place makes more neutron stars than a metal-poor one. At the middle draw a standalone sector rolls 4.08e-3 per star (Boss's research density is 4.0e-3, now a test target); bound plus ejected neutron stars equal what the galaxy formed.
- The neutron star unbound share (0.2 to 0.4) is drawn once with the galaxy's structure and read by both the bound count and the ejected field. `galaxy/black_hole_density.py` is now `galaxy/remnant_density.py`, taking the remnant kind.
- Black hole counts are unchanged. `docs/design/neutron-star-density.md` records the model.

generation-output: changed

TODO items resolved: GEN.554 slice 1 (neutron star counts follow the mass function; GEN.554 stays open for the mass islands and the kick step) / Added: none
