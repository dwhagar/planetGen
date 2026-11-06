**Headline:** The Sector Map draws every nebula and supernova remnant as a soft volume of points, with no hard ball or ring left.

### Changed
- **Supernova remnants and nebulae with no stored id are drawn from points on the Sector Map.** A remnant used to show as a hard-edged analytic ball and shell, and a scattered nebula or one a waiting star lights as a smooth ball, because only a stored nebula could ask for its points by id (MAP.142). They now ask for the points of their own sphere in the colour of their family, the same sampler and density the stored nebulae use, so they match the Galaxy Map (MAP.392). The analytic sphere stays only as the fallback of a host that cannot draw points.
- New route `GET /api/galaxy/cloud-volume?x=&y=&z=&radius_pc=&family=` (and the site's `/galaxy/cloud-volume`) returns the packed points of a sphere; additive, no existing route changes. `generation-output: unchanged`.
