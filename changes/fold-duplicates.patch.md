**Headline:** Twelve pieces of code that were written twice are written once: the body child tables, the block chain, the stored-tree read, the stars filter check, the admin check, the HTML and JSON escapes of the 3D map, and the rest of the list the layered-fixes read left (GEN.1691, UX.312). Pages and generated output are unchanged.

### Fixed
- **A bad tier or tweak code on the References pages is an error, not a blank (GEN.1691).** `references_site._lookup` raises `KeyError` for a code the legend has no row for, where it returned an empty name and the page showed a blank.

### Changed
- **Storage and generation (GEN.1691).** `fingerprint.BODY_CHILD_SUFFIXES` is the one list of the child tables every body kind has, and an edit deletes those rows by it (`db/edits`). `galaxy/drill.chain_blocks` and `STORED_LEVELS` are the one block chain `db/block_stats` and `db/block_lights` walk; `block_stats.clear` empties its table in batches as `block_lights.clear` does. `scatter_systems.read_member_tree` reads through `read_stored_trees` (whose `StoredMember` now carries a group member's `offset_au` and `velocity_kms`). `detmath.float_or` is the NULL-to-float helper of `db/sky_sphere` and `db/sky_sprites`. The belt and field member pages share `belt_member_info.neighbour_designations`, and the field pages build their field block once (`_field_block`).
- **Web (UX.312).** `db/route_neighbours` reads the nearest ladder from `db/near` (`ladder_pc`) instead of its own copy. The Generate pages, the job page and the system page check the admin with `admin_pages._require_admin`; `generate_page._admin_or_redirect` is gone. The center fact of a nebula or remnant class is a cloud fact (`classref._cloud_facts`). The life-air chart picks each gas's series class once. `fmt.json_script` is the one `<script type="application/json">` escape (the 3D map and the JSON-LD of the front page), and the 3D map escapes text with `fmt.esc`. `util/star_count.valid_text` is the one stars-filter check of the Nearby and Search pages.

generation-output: unchanged
generation-affects: none
