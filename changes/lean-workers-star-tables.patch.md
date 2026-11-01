**Headline:** a worker holds the star sampling tables only for the stages that draw stars

### Changed

A worker used to build the sampling tables of every star population when it started, about 1.4 GB, so the stages that draw no star (the link, the settle, the map-luminosity fill, the habitability re-score) ran in workers of 1.5 GB each. Now the units that draw stars build the tables in the worker before their horse forks, and a worker that holds them starts itself over, putting the unit back at the front of its queue, before a unit that needs none. A worker and its horse measured 601 MB for a link (1,564 MB before), 588 for a settle, 552 for a re-score and 540 for a map-luminosity fill, and the link ran in 57 s instead of 138 s. The memory figures that hold a stage to what the host's memory holds are now 700 MB for these stages (11 workers on a 16 GB host, 4 before), 1,700 MB for the scatters and 2,600 MB for the sector stages. No generated output changes.
