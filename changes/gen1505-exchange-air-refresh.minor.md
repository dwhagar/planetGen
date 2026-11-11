**Headline:** A planet whose eccentricity the secular exchange changes, however little, now has its air and temperature made again on the new orbit, giants included.

### Changed

- After the secular exchange of GEN.1450 moved a planet's eccentricity, its air, class and temperature were made again only for a change over 1e-4 (`ECCENTRICITY_SETTLE_TOLERANCE`), and a giant's only if its radius moved, so a smaller change left the mean light, the temperature and the air on the old orbit (a sweep body at 310.788 K against 310.774 K). Any change of the eccentricity now marks the planet's air stale (`planets.stale_air`), and a giant's air is made again in `settle_giant` on the orbit it has then, whether or not its radius moved (volume and gravity are only remade when it did). A planet the exchange leaves alone is unchanged.

generation-output: changed
generation-affects: fill
