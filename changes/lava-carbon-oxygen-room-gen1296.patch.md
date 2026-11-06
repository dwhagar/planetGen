**Headline:** A hot lava world that holds carbon in its air no longer stops the build with a solver that does not converge.

### Changed

- The carbon gas of a volatile-bearing lava world (`physics/mixture_lava.py`, row `i6`) is solved at a gas pressure that can hold the oxygen the redox buffer demands (GEN.1296). The buffered oxygen (O2 and O) alone filled the surface pressure from about 2600 K at a millibar or 4000 K at a bar, so `gas_equilibrium.solve` had no room for carbon and raised `no convergence`: 282 of 1200 draws of a sweep from 1400 to 6000 K at 0.001 to 100 bar. Where the oxygen takes under 99 percent of the surface pressure (`OXYGEN_SHARE_MAX`, stand-in row AM96, tier 6, tweak level 0) nothing changes: all 918 draws of the sweep that converged before give the same mixture to the last digit. Where it does not, the gas pressure is what holds the oxygen with 1 percent to spare, and the carbon gases are scaled to it.

generation-output: unchanged
