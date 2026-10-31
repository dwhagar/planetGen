**Headline:** The Galaxy Map's history entries keep only the camera's address; Back, Forward and a reload work the stage out from where the camera is.

### Changed
- **History entries no longer carry a stage (MAP.337, steps 3 and 4).** `galaxystageview.js` stores `{galaxyStage, mapIndex, maxIndex}` in each entry and derives the stage the page loads on from the `cam` address (`stageForCamera`), where an entry used to hold the stage it was made at. Dead code is removed (`fovHalf`, `optionAt`, `matchAddress`).

No generator output changes.

generation-output: unchanged

TODO items resolved: MAP.337 (steps 1 to 4 merged; `parentStage` stays as the container trail the crumbs and Up read, MAP.227).
