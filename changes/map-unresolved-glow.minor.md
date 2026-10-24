**Headline:** The Galaxy Map shows the light of the stars too far away to tell apart as a soft glow that turns into stars as you zoom in, and the Menu can switch it off.

### Added
- **The unresolved glow (MAP.266).** Each sector-block node lists only its brightest stars; the rest of its light, stored in its aggregate (the pre-placed stars nothing was built from and the generated sectors' total), is drawn as one soft blob at the centre of that light, in the color of its temperature, under the stars. A blob shows for the share the cut gives its node (`galaxynodecut.weightCalculator`), so as a block hands over to its children the glow gives way to their stars and blobs and the light is conserved: nothing pops and nothing doubles. `static/galaxyglow.js` chooses the blobs and how strongly each shows (a tone curve over its apparent magnitude, judged 5 magnitudes deeper than a point of light); `static/galaxyglowlayer.js` draws them as screen-sized quads. The Menu's "Show on the map" group has an "Unresolved glow" button (on by default; `hide=glow` in the address turns it off). The Sector Map has no glow.

### Tests
- A second warm-up over addressed sectors builds no answer (MAP.289: the bug does not reproduce; the test holds it).
