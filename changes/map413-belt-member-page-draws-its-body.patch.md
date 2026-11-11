**Headline:** The page of a belt member now shows a picture of the body: a lump or a ball in the colour of its taxon and makeup.

### Changed

- MAP.413: `/object/belt:<ID>/member/<rank>` has a "Sensor view" panel above the table. `bodyview.member_view(designation, member)` makes the look record with `bodylook.member_look` (MAP.412) from the page's member block, its seed from the designation so the 3D view (MAP.414) draws the same body; the existing `static/bodyview.js` paints it (a body under the roundness limit through `lumpyview.js`, a round one as a sphere). One body per page; nothing is guessed, so there is no estimate line. The caption says the colour comes from the taxon and the shares of rock, metal, water and ice.
- Reused: `bodyview.js`, `lumpyview.js` and `lumpymesh.js` unchanged, the item page's Sensor panel markup, `bodylook.member_look`. Checked in headless Chromium: a lumpy dark carbonaceous body and a round stony one paint with no console message.

generation-affects: none
