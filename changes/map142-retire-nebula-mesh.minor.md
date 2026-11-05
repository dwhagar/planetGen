**Headline:** The hard-edged nebula shell is gone from the maps: nebulae are soft clouds everywhere.

### Changed
- The Galaxy Map, the Sector Map and the nebula page no longer draw a nebula as a translucent hard-edged mesh. A nebula is its own soft cloud of points (or, where the browser cannot draw that, the soft sprite or sphere it had before the mesh). The site no longer serves `/galaxy/nebula/<id>/shape`; the API's `GET /api/nebulae/<id>/shape` is unchanged.

generation-output: unchanged
