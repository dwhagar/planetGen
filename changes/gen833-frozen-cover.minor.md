**Headline:** A cold moon that holds its water as ice is covered by the ice it holds, not by a toy ocean cover.

### Changed
- **GEN.833: the ground draw of a frozen body reads the real ice cover.** A body in the state `ice` or `ice-covered ocean` now stores as its ocean fraction the cover of all the water it holds at ice's density (`hydrosphere.frozen_cover`), with the relief capped at 10.5 km (observed icy relief is a few to 20 km at any gravity), not the cover of the part of its water a liquid ocean would have had. A 200 km moon with a fifth of its mass as water used to draw bare basalt over two thirds of its surface under a 130 km ice shell and is now all ice; a rocky world with little water keeps a partial cover. The planets and moons whose results move are the ones in those two states below Earth's gravity or with mantle-buffered water.

generation-output: changed
