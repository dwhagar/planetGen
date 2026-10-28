**Headline:** The Galaxy Map marks every charted region as a soft sphere colored by how much of it is charted, names it, frames the region you are in, and lists them all in the Menu.

### Changed
- **Regions are what the whole-galaxy view shows (MAP.316).** Every region with at least one generated sector, however little, is a soft sphere at its extent: cool blue for a lone sector, teal, then gold and amber the more of it is charted. The spheres fade as you come in and the region's own stars take over. No sector box is drawn with them. A layer of the Menu ("Charted regions"); `hide=regions` in the address turns it off.
- **Clicks on a region.** Hovering a sphere names it, a click shows how much of it is charted with a Fly to region button, and a double click flies to it.
- **A charted-regions list in the Menu.** One button a region ("Fardoo Beez Region, 14 percent charted"), sorted by share charted or by name, the always-available way to reach known space.
- **A frame round the region you are in (MAP.327).** In the galactic view the same locked wireframe that frames an opened sector stands round the charted region the view is in, with nothing dimmed.
- **New endpoint.** `GET /api/galaxy/charted-regions` lists the regions with their extents, shares and names.
