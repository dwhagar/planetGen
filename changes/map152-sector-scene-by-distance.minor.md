**Headline:** On the Galaxy Map a sector's own stars and phenomena fade in as the camera comes up to it, with no click and no change of view.

### Changed
- In a window of sectors, a charted sector's scene (every star, nebula and body of it) is fetched once the camera has stayed within 1.5 sector edges of its centre for a quarter of a second, fades in as the camera comes nearer, is whole at 1 edge, and is dropped past 3 edges. The sector's block fades out as its scene comes in. Opening the sector (a click, or one more wheel tick) takes the scene over as it is: no second fetch.
- A window of sectors lets the camera come to half a sector edge from the sector at its middle (it stopped about 1.4 edges out, short of where the scene is whole). One more wheel tick at that distance opens the sector where the camera is, without a flight, and the next goes on in toward a star.
- A star of the scene that is up can be hovered, picked and opened ("Open system here") once about a third of it shows.
