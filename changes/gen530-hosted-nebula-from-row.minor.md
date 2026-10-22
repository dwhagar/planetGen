**Headline:** A pre-placed star's nebula shows before its sector is built: the star's page names it and the uncharted sector's map draws it.

### Changed
- **Hosted nebula from the star's row (GEN.530).** An O, B, A or Wolf-Rayet star that the scatter placed but nothing has built now shows the nebula it lights, the same one its sector builds. The uncharted star's page has a "Lights a nebula" line with its name, kind, class and radius; the map of an uncharted sector draws the nebula at the star, and clicking it opens the star's page. A star already built into a system shows none there, since the system holds it. `GET /api/uncharted-systems/<id>` gains `hosted_nebula` and `GET /api/galaxy/uncharted` gains `hosted_nebulae`.
- The Galaxy Map's tiles are unchanged: they list scattered phenomena by their own rows, and a hosted nebula has none until it is built.

TODO items resolved / Added: GEN.530 / Added: none
