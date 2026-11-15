**Headline:** The cosmogenic isotopes of a planet's surface (10Be, 26Al) are now made with the planet's magnetic field, which the first record was made without.

### Changed
- GEN.1544: the isotope record is first made with the air (`planets.calculate_atmospheric_conditions`), before `magnetism.generate_field` draws the dipole, but the surface cosmogenic rows (GEN.1489) are cut by the dipole moment. A planet or moon with a field stored rows made as if it had none until the final-age pass ran them again. `generate_spin` now makes the record again right after the field (own streams, no draw of the run's moves). Planets with a field store smaller 10Be and 26Al surface amounts, by up to the field factor (about 0.59 for an Earth-size dipole). `tests/test_isotope_record.py` pins a planet with a field: its stored record is the one the finished body gives.

generation-output: changed
generation-affects: fill
TODO items resolved: GEN.1544
