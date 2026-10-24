**Headline:** Rogue planets, comets, black holes, neutron stars and asteroid fields are now drawn from their own stored data instead of one fixed look each.

### Changed
- **Rogue planets use the planet painters (MAP.278).** A rogue planet's View is painted from a look record built from its stored row (`bodylook.rogue_planet_look`): a rocky one by the solid-world painter (ice where its surface is cold, ice over an ocean when it is stored as liquid under an ice shell, no sea otherwise), a giant as bands of its temperature colour. It replaces the fixed blue-grey.
- **Comet nuclei take a colour from activity.** An active comet's nucleus is grey-white, an inactive one dark red-brown, and its lumpy shape is the same on every visit.
- **Black hole disks and neutron star surfaces take their temperature colour.** The accretion disk's hot colour and the neutron star's surface are the blackbody colour of the stored temperature, in place of one blue-white or yellow-white for each class.
- **Asteroid field grains follow the composition.** On the Sector Map an asteroid field's cloud and glow are tinted by its composition family (carbonaceous dark, stony tan, metallic grey, icy pale, basaltic brown), from the field's class letter.
- **Each phenomenon View has its own stars.** The background starfield is seeded from the object, not one fixed seed for every View.

Nebula tints by type already existed; supernova remnant shells by age are not in this change (they need the remnant's age on the sector scene).
