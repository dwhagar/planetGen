**Headline:** Tests: a double click on a star or a body flies to it and opens its information

### Fixed

- MAP.140 (GitHub issue 714): the behaviour was built by the fly-through (MAP.146, MAP.226, MAP.304) and the select-and-open that goes with it was only half covered. Two browser tests now hold it: on the Galaxy Map a single click on a star selects it and opens its information without a flight, and a double click flies the camera onto it with the information still open; on the System Map's 3D view a double click on a planet flies the camera closer to it and its name is in the information panel.
