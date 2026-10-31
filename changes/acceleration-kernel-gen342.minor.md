**Headline:** The pull on an object from the masses around it and from the galaxy is now one function, with the velocity kick that uses it.

### Changed
- **The acceleration kernel (GEN.342).** `planetgen/galaxy/acceleration.py` returns the acceleration at a point as the pull of the object's influence set (the exact nearest masses and the folded monopoles, unsoftened but for 0.001 pc that keeps two masses at one point finite, with no minimum-acceleration cut) plus the galactic gradient of the potential, and turns it into a new velocity (`kick`, `velocity_update`). A star takes the whole pull; a planet or moon takes the pull at itself less the pull at its host (`host_frame_pull`), so only the perturbation is left. It is a pure function, so the same masses always give the same pull. Nothing calls it yet: the star stepper (GEN.340) and the orbital update run (GEN.105) do, so generation output does not change. Also corrected the time unit in the potential's orbit test: a Myr is 1.0227 pc per (km/s), so the five-gigayear orbit it checks now is five gigayears.

generation-output: unchanged
