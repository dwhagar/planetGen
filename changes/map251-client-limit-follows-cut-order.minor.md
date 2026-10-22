**Headline:** When the Galaxy Map has more stars on screen than it can show, it dims the satellite-galaxy stars first, then the halo's, then the farthest and faintest.

### Changed
- **The on-screen star limit follows the cut order (MAP.251).** The limit that holds about 20,000 stars on screen (8,000 on a phone) used to judge every star by its light alone. A star of a satellite galaxy now counts as 6 magnitudes dimmer and a halo star as 3, both when the page counts the stars on screen and in the shader, so a crowded view dims those first, then the rest by apparent magnitude: the nearest and brightest go last, and the farther away an object is the sooner it goes. A view with room for every star shows them all.
- The stage views and the Sector Map spend no star budget (a stage draws stored block counts, the Sector Map every system of one sector), so the rule is written down for them and nothing changes there.

TODO items resolved / Added: MAP.251 / Added: none
