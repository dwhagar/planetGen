**Headline:** A planet draw no longer hangs or crashes on a world whose surface is hotter than its interior would be.

### Fixed
- **A planet draw around an F5V star at an ecosphere orbit never returned (GEN.866).** The thermal history of a wet rocky planet with a 2300 K surface let a step overshoot below zero kelvin; the viscosity law on that temperature gave a flux of 1e59 W/m2 and a step of 1e-65 Gy, so the run never ended. The mantle temperature is now held between the surface temperature and the ceiling after every step, as the rate already was. No value on a normal path changes.
- **A system draw raised ZeroDivisionError in the interior heat for about 2 in 600 F5V and K0III systems (GEN.895).** The same overshoot gave a mantle viscosity of exactly 0; the same clamp removes it. The census repro now draws 300 F5V systems without the error.
