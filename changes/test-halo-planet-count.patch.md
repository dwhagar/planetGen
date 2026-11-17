**Headline:** The halo planet-count test counts the planets metallicity scales, not the class D tiny bodies it never touched.

### Fixed
- The test that a metal-poor halo star carries fewer planets than a disk one failed on main (halo mean 17.9 against disk 17.4 over the same 30 systems each) because class D tiny bodies, 13.2 a system in the halo and 9.8 in the disk, were counted in `len(system.planets)`. Without class D and belts the means are 3.0 and 6.1, a ratio of 0.50, and the 0.6 bound is unchanged.

generation-output: unchanged
