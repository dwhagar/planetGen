**Headline:** The lists, the counts, the search and Nearby now include the scattered black holes, neutron stars, remnants, planetary nebulae and hypervelocity stars and the field clouds that no sector has stored yet, so the Lists page counts what the map draws.

### Fixed

- The object index reads only the tables that hold a stored object, so the Lists page, the admin Object index table, the Phenomenon facet and panel of the search page and Nearby showed NEBULA 0, SUPERNOVA_REMNANT 0 and BLACK_HOLE 1 while the map draws the 110 thousand field clouds (`nebula_field`) and the 9.7 million scattered phenomena (`phenomenon_scatter`). They are not put in the index (about 3 GB, three times the cache budget the index is kept first in). `db/uncharted.py` reads them beside it: a name-sorted list shows the charted objects first and the uncharted ones after them (field clouds, then scattered rows), and a count adds the two.
- A scattered row links to its scatter page (`/uncharted-phenomenon/<id>`); a field cloud, which has no Object ID or page until a sector stores it, is listed by its key (`field-i.j.k.n`).

### Added

- A Hypervelocity stars list (`/lists/hypervelocity_star`, `GET /api/object-lists/hypervelocity_star`), the Hypervelocity Star type in the search page's Phenomenon facet, and the kind in Nearby: a hypervelocity star is a system once its sector is built, so until then only the scatter row has it.
- Nearby lists the uncharted objects in range (`uncharted`, `scatter_id`, `ref` null) instead of counting the scattered ones under `not_generated`; the pre-placed bright stars are still counted only.
- The admin Object index table has an Uncharted column beside Rows, recorded by each index run.

generation-output: unchanged
