**Headline:** An object's page shows its Sensor view first, with the "Sensors show ..." line under the picture.

### Changed
- **The Sensor view panel of the object page** (UX.276): the picture of a star, planet or moon is now chosen through the Sensor registry (`web/lib/sensor.py`, UX.275), and the "Sensors show ..." line of the description sits under it instead of in the notes further down. Where the picture guessed an input (a world's ground color or a giant's cloud-top color, worked out from the temperature because no ground mixture or cloud chemistry is stored) the panel says it is an estimate and why.
- The picture on a phenomenon's page (a neutron star, black hole, quasar, rogue planet, interstellar comet or nebula) is headed "Sensor view" instead of "View". A belt or a comet of a system has no picture yet, and so no panel. `generation-output: unchanged`.
