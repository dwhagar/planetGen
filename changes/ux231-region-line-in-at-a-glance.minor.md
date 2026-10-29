**Headline:** Sector, system and object pages now say which charted region they lie in, such as "in Vang Zarh Region".
### Added
- **Region line (UX.231).** The At a glance block of a sector, a star system and every star, planet, moon, belt, comet and phenomenon page has a "Region" row, "in Vang Zarh Region", when the position lies in a charted region (a block of nine that holds a generated sector). A position in no charted region, or a galaxy with no naming key, shows no row. The API answers for a sector, a system and an object carry `region` (`name` and `parent`) for it.
