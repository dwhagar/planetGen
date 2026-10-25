**Headline:** The Galaxy Map names the stars and other objects it draws, and a click on a name selects that item and flies to it.

### Changed
- **Names on every visible item (MAP.298).** Beside the sectors, systems, planets and nebulae it already named, the Galaxy Map now writes the name of each named star it draws and of each named object of the other kinds (comets, rogue planets, black holes, neutron stars) where there is room, the brightest first. Two new entries in the Menu's Labels group, "Star names" and "Object names", turn them off.
- **A name is a link.** The pointer over a name marks it; a click on it selects the item (its details open in the map's window) and flies the camera to it.
- **Packed nodes carry star names.** The packed node (`TILE_FORMAT` 4) lists the names of its named stars by id, so the cached nodes are rebuilt once.
