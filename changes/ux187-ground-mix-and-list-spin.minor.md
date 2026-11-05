**Headline:** The page of a planet or moon now shows what its ground is made of, and the lists show each body's spin, tilt and internal heat.

### Changed

- The Composition block of a planet or moon page lists the ground type that covers most of it, the grain of the ground, and the share of each of the 13 ground types (ocean, ice, ultramafic, anorthositic, basalt, granitoid, intermediate, clay, iron-oxidized ground, metal, carbon, evaporite salts, sulfur), each with its unit and a link to the method page that gives the source. The property sheet of the API (`GET /api/objects/<ref>/sheet`) carries the same rows.
- The planet and moon lists have three more columns, rotation period, axial tilt and internal heat flow, and the star list has the rotation period; the object-list API (`GET /api/object-lists/<kind>`) returns the same values on each item.
