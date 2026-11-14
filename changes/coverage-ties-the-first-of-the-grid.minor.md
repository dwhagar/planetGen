**Headline:** The hottest and coldest place of a world and its best site no longer flip hemisphere on the last bits of a float.

### Changed
- **Ties in the temperature coverage (found by the habitability re-score test of Train 26).** A world whose climate is the same north and south (an eccentric orbit with the spin axis in the plane of the perihelion) has its hottest sample, and its best site, at two mirror places that agree to the last float bits. `temperature_coverage._extreme` took the larger by an exact `argmax`, and the best site broke its tie on the exact float swing, so the stored world and its re-score from the database could name opposite hemispheres. The hottest and coldest are now the first sample in the grid's order within a millionth of a kelvin of the extreme, and the best-site tie goes by the swing rounded to a millionth, then the larger blue-band share of the time, then the mean nearest the middle of the blue band, then the grid's order. A stored best site on a world that had such a tie can name another cell.

generation-output: changed
generation-affects: fill
