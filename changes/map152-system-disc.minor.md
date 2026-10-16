**Headline:** On the Galaxy Map a star's system fades in as the camera comes up to it, and opens from that scene with one more wheel tick.

### Changed
- In an opened sector, a star's system (its planets, moons and orbits) is fetched once its disc would be about 20 pixels across, fades in to 40 pixels, and is dropped below 10 pixels. It is not opened by that: the address and the panel stay the sector's, and one more wheel tick at the nearest zoom opens the system from the scene already up, with no second fetch.
- In an opened sector the camera can now come right up to a star near the middle of the view (to where its disc is 60 pixels across), instead of stopping at an eighth of the sector's fit. Move the middle onto the star first (Shift-click or right-click on it).
