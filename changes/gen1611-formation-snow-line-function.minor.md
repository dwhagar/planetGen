**Headline:** The snow line a star's disk was built round can now be worked out for an evolved star from its birth, not from the luminosity it has today; nothing reads it yet.

### Changed

- GEN.1611 slice 0: new `physics/formation_snow_line.py` (`formation_luminosity_sol`, `formation_snow_line_au`): a main-sequence or contracting star keeps today's luminosity; an evolved one gets its track's luminosity at `FORMATION_DISK_AGE_GY` (0.01, stand-in PA18). A white dwarf born at 4.36 solar masses gets a line of about 42 AU where today's luminosity gave 0.02 AU. The generator does not read it yet (wiring register); slice 1 wires it, with the sourced scatter cap, as a lock item.

generation-output: unchanged
