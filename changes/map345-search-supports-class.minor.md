**Headline:** The search page can now find planets and moons that could support a given kind of life, with a minimum score.

### Changed
- **Search by the life classes a world supports (MAP.345).** Two new filters on the search page, one for planets and one for moons, pick the worlds whose conditions support a given life class (1a, 1b, 2, 3a and so on) at a score of 1 or more or 5 or more out of 9. They sit beside the filter for the class that actually lives on a world, and each option shows how many worlds it holds. Worlds with no class scores (gas giants and worlds made before the scores) match none. The same filters work in the API (`supports=1a:5`, `moon_supports=2:1`).
