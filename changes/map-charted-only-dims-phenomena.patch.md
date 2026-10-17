**Headline:** "Charted only" on the Galaxy Map now dims the black holes, neutron stars and quasars of sectors that are not filled yet.

### Fixed
- With "Charted only" on, the scattered black holes, neutron stars and quasars outside the filled sectors stayed fully lit while the stars dimmed, because the map treated every phenomenon as being in a filled sector. A scattered point now carries `charted` (true once its sector is built) and the map dims the unbuilt ones the way it dims stars outside charted sectors; phenomena of built sectors stay lit.

TODO items resolved: MAP.252 / Added: none
