**Headline:** A planet's or moon's Life section now has a Life and air chart: how its oxygen, carbon dioxide and methane changed as the life took hold.

### Changed
- **The Life and air chart.** Under the Life section of a planet or moon whose life changed its air there is a collapsed block, "Life and air through the ages". It plots the partial pressure of oxygen, carbon dioxide and methane against the age of the star, with the stages of the life (first life, producers, complex cells, animals) numbered on the chart and the same figures in a table. The air is worked out backwards from the air the world has now, when you open the page; nothing new is stored. A class with tuned numbers (Classes 6 to 9 and others at tweak level L3 or L4) is marked "Speculative chemistry" with its tier. A world whose life did not change its air shows no chart.
- **`GET /api/planets/<id>/life-air` and `GET /api/moons/<id>/life-air`** give the chart's data (new routes; the API version does not change).
