**Headline:** The Life and air chart and its route now sum every life class that lives on a world, and the page says which class made how much of each gas.

### Changed
- **GEN.1213 item 3: `GET /api/planets/<id>/life-air` (and moons) drew the history of the dominant class only.** It now reads the classes stored for the body (`planet_life`, `moon_life`) and sums them in timeline order (`life_air.air_history`), so `classes` lists every class with its share of the food production. Each gas also carries `by_class`: the share of what life changed in that gas that each class made (the methane, and the carbon dioxide it is made from, by each class's methane flux times its share; the oxygen by the oxygenic classes' shares).
- The chart's section of the planet or moon page lists those shares under the legend when more than one class changed a gas ("Methane: Class 3 75 percent, Class 1a 25 percent"). A world with one class reads as before.
- No generator or stored-data change: the history is worked out on request.

TODO items advanced: GEN.1213.
