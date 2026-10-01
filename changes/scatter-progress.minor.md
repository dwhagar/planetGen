### Changed

- The bright-star scatter's progress bar and its time left now count each layer's expected stars (a quick sample of the density model before the scatter starts, plus a little for every ring walked) rather than layers, and move as layers in progress report their stars, so the thin layers at the edges of the disk no longer throw the estimate off (PERF.9). It shows a share done, "Bright stars (12 of 1,271 layers) 34%".
- While layers take longer than 30 seconds each, a second bar under it shows the stars of the layers being drawn, done of their estimate, with its own time left; it goes once layers finish faster than one every 20 seconds (PERF.4). The Generate page shows both, and each layer's time goes into the speed stats as a "scatter" task.
