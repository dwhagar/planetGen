**Headline:** Code that read or picked by a planet class letter before the class is decided is removed; the Hycean box and the trojan room test now read mass, radius and temperature alone.

### Changed
- **Class-letter leftovers removed (GEN.1400).** `trojan_bodies.classes_allowed` and its two class constants (rows TB1 and TB2) are gone; `fits` and `room_earth` test the mass alone. `hycean_candidate` takes mass, radius and equilibrium temperature and no longer needs the letters T, I or V. The look of a body no longer carries `cls`, the air colour is not handed the class, the habitability world no longer carries `planet_class`, and the orphaned docstring about a removed milestone cap map is deleted. `generation/system.py` `m_count` waits for GEN.688.

TODO items resolved: GEN.1400 (all but `m_count`)
