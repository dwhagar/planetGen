**Headline:** A young star can now be given its own disk dispersal time, drawn once from the survey curve, ready for the build that keeps its disk until then and its planets until after.

### Changed

- GEN.1551 (built ahead for Lane 11): new `physics/disk_dispersal.py` with `disk_dispersal_time_yr` (a star's own time, from a stream keyed on its mass and a birth key and never its age; 0.8 of disks left at 1 Myr, 0.5 at 3, 0.1 at 5.5, none by 10; a companion inside 47 au cuts two thirds of disks inside 1 Myr), `disk_present`, `disk_share` and `survival_to_time_yr`; stand-ins DD1 to DD3 (registry 6.274). Nothing calls them yet; GEN.1551 wires the star's disk stage and GEN.1552 the planet cut.

generation-output: unchanged
generation-affects: none
