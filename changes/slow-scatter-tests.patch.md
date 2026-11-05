**Headline:** Three bright-star scatter tests that took four to nine minutes under load now take one to two; no generated output changes.

### Changed
- **Scatter tests drawn only as far as their checks need (TEST.417).** In `test_bright_star_scatter.py`, the test that checks backfilled cells hold band stars in their own cell now takes three cells of the dense shell 2 (about 120 stars each; twelve took nine minutes) and keeps the twelve of shell 3; both shells are still checked for the cell, the band and the Poisson count. The bar test for the "Topping up backfilled sectors" phase now takes the band down to 400 solar luminosities (the cheap first step, which reports the same phase) instead of down to the floor.

generation-output: unchanged

TODO items resolved: none (TEST.417 continues).
