**Headline:** the nine uses of `math.hypot`, `math.dist` and `@` the reproducibility scan found now go through `detmath`, so a galaxy's bits no longer depend on the Python or BLAS library that made it.

### Fixed
- `test_no_version_or_library_dependent_maths_in_the_package` was red on main: `generation/system.py` (`math.dist` of the dark members' places), `physics/eclipses.py`, `physics/influx.py` and `physics/kozai_planets.py` (`math.hypot`), and `physics/obliquity_variation.py` and `physics/secular_chain.py` (`@` on the complex mode matrices, which runs through BLAS) now use `detmath.dist`, `detmath.hypot` and `detmath.dot`.
- `detmath.matmul` and `detmath.dot` take complex arrays. numpy's complex product is not a fixed operation (its vector loops fuse the multiplications and additions on some shapes and processors), so each product is formed from real ones, the real part `ar*br - ai*bi` and the imaginary part `ar*bi + ai*br`, and the products are added left to right. `dot` of two vectors returns a Python number (float, or complex for complex input).
- Stored values move within tolerance: a census of 60 seeded systems (782 bodies) changed in 265 bodies, almost all by under 1e-7 relative (seasonal temperature extremes, which the solver stops at a tolerance; the largest is 1.1e-7), and in one body a threshold decided differently (its ocean pH, flavor text count and water activity).

generation-output: changed within the solvers' tolerance (no golden, no epoch, lock item)
