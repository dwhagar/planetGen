**Headline:** Stars and black holes, neutron stars and quasars are points of light on the maps until you fly near a sector.

### Changed
- On the Galaxy Map and in an opened sector, a star is now a point of light about 1.5 pixels wide that does not grow as the camera approaches; its brightness tells stars apart. About one sector away (16 parsecs, one named constant shared by both scenes) it grows smoothly into its core and glow over the last two zoom steps, and shrinks back to a point flying out, so the hand-off from the Galaxy Map into a sector shows no jump in size.
- A black hole, neutron star or quasar is drawn the same way and as bright as a star of its mass: a 10 solar mass black hole looks about as bright as a 10 solar mass star, a black hole a little dimmer (0.8 of it) and keeping its purple, a neutron star keeping its dark blue. One with no mass on record is drawn as a one solar mass star, and none is drawn brighter than the brightest star on the map. The Galaxy Map's tiles and the sector scene both take their brightness from `phenomenon_map_luminosity_sol`, and the tiles now carry each point's `mass_solar` and `map_luminosity_sol`.
- The system map is unchanged: it has no background stars of its own and draws its bodies at their drawn size.

TODO items resolved: MAP.208 / Added: none
