**Headline:** Three tests that were red on main are repaired: the sky star rows have their burn-down lines, the Large Cloud's supernova claim follows the new solar circle, and the halo planet count is measured against the disk's.

### Fixed
- `test_burn_down_register`: the nine hand-tuned sky rows of section 6.181 (SKL1 to SKL9, the look of a star, the star caps and the page's reach) have register lines saying what would replace each.
- `test_gen_magellanic_classes` (found by bisecting to GEN.1230): the Large Cloud's supernova remnant factor was 1.23 and is 0.92. GEN.1230 gave the thin disc an inside-out infall, which raised the solar circle's young share from 0.011 to 0.016 while the Cloud's stayed 0.014 (a 0.39 solar masses a year Cloud, about what is observed). The claim "above 1" belonged to the old reference; the test now checks the Large Cloud's factor is near 1 (between 0.5 and 2) and larger than the Small Cloud's, and the cell still equals the component pass's factor.
- `test_gen_kind_factors` halo planets (found by bisecting to GEN.445, epoch 74): a disk G dwarf holds about 9.3 planets now, not the typed 6.9, and a count of planets is heavy tailed (a 17-planet halo system sank the 20-system mean). The test measures the disk on the same draws and requires the halo's mean (about 3.4 to 4.3) to be under 0.6 of it.

### Notes
- The nebula group of TEST.580 stays with Tides 1. Train 10 changes `test_gen_kind_factors.py` in other tests (rogue planets), not this one.
