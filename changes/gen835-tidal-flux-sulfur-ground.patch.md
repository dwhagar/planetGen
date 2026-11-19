**Headline:** A moon or planet the tide heats harder than its own heat now gets the sulfur ground of an Io, which no body had before.

### Fixed

- **The ground draw never read the tidal heat flux.** `surface_composition.inputs_of` did not pass `tidal_flux_w_m2` on, so the sulfur rule ("sulfur where the tide beats the body's own heat") always read zero. On 200 systems (3,256 rocky bodies) sulfur ground goes from 0 bodies to 41: 8 of 224 moons and 33 of 3,032 planets (81 bodies have a tide above their own heat; the rest are ice ground or have none, where the rule does not apply). One body changes ground type.

generation-output: changed
generation-affects: fill
