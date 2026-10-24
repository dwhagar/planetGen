**Headline:** The system page and the sector page list their contents as trees of links to each object's own page.

### Changed

- **The system page (UX.173).** The Stars table and the dropdown rows of stars, planets, moons, belts and comets are replaced by one tree: the system's stars, each already open onto its planets, belts and comets, every row a link to that object's own page, with moons and deeper levels loaded when you open a row. The text that used to drop down under a row is now on the object's page. A planet's or moon's habitability chips, a life world's dominant species and the facilities a body hosts moved to that body's page, and a body's Admin menu is there too (the system page keeps its own).
- **The sector page (UX.174).** The arrow before a system's name in the Contents table opens its stars and bodies under it, as the same tree; the table's sorting, filters and paging are unchanged.
- **Links to a body (UX.176).** Search results, descriptions and map popups all open the body's own page; the old `/system/<id>#planet-12` row links, and the script and styles that opened and highlighted those rows, are gone.
