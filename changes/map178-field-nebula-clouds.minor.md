**Headline:** The Galaxy Map now draws the galaxy's dark molecular clouds one by one, before any sector near them has been generated.

### Changed

- A map block of 81 sectors or fewer lists the clouds the galaxy's seed places (the `nebula_field` table, classes M and N, giant molecular and dark clouds) that no sector has stored yet, beside the stored nebulae, and the map draws each as a soft dark-dust sprite that fades in between 1 and 4 pixels of radius (MAP.178). A cloud a sector stores is drawn once, as the stored nebula. A block larger than that keeps showing the nebula regions (MAP.176), which hold these clouds.
- A field cloud has no page, name or stored shape yet, so it is not labelled, picked or opened as a mesh; picking and fly-to for clouds and regions is MAP.180. The map's saved blocks are rebuilt once on the next update.
