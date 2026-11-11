**Headline:** The 3D system view now draws the heaviest bodies of each asteroid belt one by one, as lumps or balls in their own colours, among the belt's dust.

### Changed

- MAP.414: `GET /api/systems/<id>/scene` gives each belt a `members` list: the heaviest `BELT_MEMBERS_DRAWN` (200) bodies it stores (`systemscene._belt_members`), each with its designation, rank, radius, `look` (`bodylook.member_look`, its seed the designation's, so the page of the body (MAP.413) and this view draw one and the same body) and `position_km`, where its stored orbit puts it at the epoch (`system_bodies.relative_position_au`). A body that has since become a moon is left out.
- `static/systemview3d.js` hangs each from what the belt goes round: a lump (MAP.279) when it is under the roundness limit, a sphere otherwise, painted by the same texture tiers, queue and cache as the planets. They are not pickable, labelled or dotted (a belt keeps hundreds); the dust points stay. `static/systemscale.js` gains `layout.radiusFor(kind, km)` and a `belt_member` radius: smaller than any moon when compressed, true size otherwise.
- Reused: `member_look`, `lumpyGeometry`, `paintSolid`, the texture queue and cache, `relative_position_au`. Tests: `tests/js/systembeltmembers.test.mjs`, `tests/test_system_scene.py` (the query and positions without a database, and with one).

generation-affects: none
