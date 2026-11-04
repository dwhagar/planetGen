**Headline:** A single nebula's cloud of points can now be fetched on its own, ready for the Sector Map and the nebula page to draw the same irregular shape the Galaxy Map shows.

### Changed

- New `GET /api/nebulae/<id>/volume` (and the web `GET /galaxy/nebula/<id>/volume`) sample one nebula's cloud from the model's gas density with the same sampler as the Galaxy Map's node volumes. Nothing is stored and no map draws it yet (MAP.142 part 1).

generation-output: unchanged
