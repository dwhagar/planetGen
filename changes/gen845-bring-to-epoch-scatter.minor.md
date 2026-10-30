**Headline:** Bringing a galaxy to the running generator epoch now draws a scatter that is behind again, instead of stopping and telling you to do it by hand.

### Changed

- **A stale scatter is redrawn by the bring-to-epoch run (GEN.845).** When an epoch has changed a scatter (the massive stars, the bright stars, the phenomena or the nebula field), `planetgen bring-to-epoch` and the Generate page's task no longer refuse. The run empties the sectors it will generate again (a scatter leaves a filled sector out), clears the stale scatters and draws them with the settings they were drawn with (luminosity floor, stellar mass limit, compact object limit), the massive stars before the bright stars that are split from them, then generates the emptied sectors again. A sector left alone because of an admin edit keeps the objects it was made with. If a scatter fails the emptied sectors stay empty and the run says how to finish by hand. Nothing a generator draws changes, so there is no epoch.
