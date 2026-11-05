**Headline:** The Galaxy Map's Color by menu has a new choice, Nebula cover, that shades each block by the share of it inside a nebula.

### Changed

- Color by gains "Nebula cover" (MAP.177). Each block or window cell is drawn in one of four steps of a single violet hue that lightens with the share of the cell inside a nebula (under 0.5%, 0.5% to 1.5%, 1.5% to 3%, 3% or more, the bins of the nebula density research), whether or not the cell holds a charted sector. A cell with no nebula is drawn as plain unfilled space.
- The steps are fixed shares, not ranked against the cells in view, so the same share has the same colour on every screen. The legend names them and the hover text of a block says its share ("0.60% inside a nebula"), or that none is charted for a sector, which has no stored cover.
- The choice rides in the address as `color=cover`, like the other fills.
