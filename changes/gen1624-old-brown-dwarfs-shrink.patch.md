**Headline:** An old brown dwarf of the lightest masses is about one Jupiter radius across, not the 1.5 to 1.7 it kept from its youth, and it is no longer up to 2.4 times too bright.

### Fixed

- GEN.1624: the BHAC15 table loses the lightest brown dwarfs while they are still young and large (13 Jupiter masses after 0.01 Gy), and `baraffe_isochrones.cooling_track` held their last radius for ever, so a 13 Jupiter mass body stayed at 1.7 Jupiter radii at 10 Gy and its luminosity (`R^2 Teff^4`) ran up to 2.4 times high. Past the age the table loses a mass its radius now shrinks as the radius of the lightest mass the table still holds does (an old brown dwarf's radius is nearly the same at every mass, 0.09 solar radii), joined to the next heavier mass's own change as before so the track still has no jump in mass (GEN.1137). A 13 Jupiter mass body at 10 Gy is 0.91 Jupiter radii. Stand-in LM4 now says so.

generation-output: changed
generation-affects: scatter, fill
