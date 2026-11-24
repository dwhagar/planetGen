**Headline:** Zooming out of an opened sector with the wheel now carries on at the same steady step instead of snapping the camera out to the whole galaxy.

### Fixed

- Leaving an opened sector by the wheel took the camera from a few parsecs to the whole-galaxy distance in one tick (MAP.257): the sector branch of the stage view's `up()` kept the camera but not its fit, so the next frame re-fitted it to the whole stage. It now keeps the fit as the system branch does, and the zoom floor no longer pushes a camera that left a system nearer than the sector's floor back up to it in one tick. A scroll from the whole galaxy to a sector, to a star's system and out again, or to an uncharted cell and out again, moves the camera the wheel's own step per tick all the way.

generation-output: unchanged
