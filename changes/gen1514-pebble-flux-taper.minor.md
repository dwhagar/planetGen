**Headline:** The pebble flux that drifts through a disk now fades as the dust-growth front moves out through the edge of the gas disk, instead of stopping dead at the edge.

### Changed

- `disk_drift` (GEN.1514, stand-in DE3): the pebble flux is multiplied by `exp(-front / edge)`, so a species is carried by a fading flux until the disk's life ends (3 Myr) and `pebble_end_myr` is gone. The core-building pebble mass has a closed form in the tapered flux.
- In the stored gas of giants, most of those of dwarf and Sun-like stars now gather their gas as the disk goes (the tapered pebbles inside the edge cannot build a 15 Earth mass core, 48 against 148 Earth masses for the Sun) and so hold few of the drifted ices; the carbon over oxygen of their envelopes is unchanged at the median and the argon, krypton and xenon fall to a few percent of the star's. Registered in GEN.725 (DE3, DE11, DE12) with the miss against the Galileo noble-gas values named as the falsifier.

generation-output: changed
generation-affects: fill
