**Headline:** The bodies an interstellar asteroid field keeps are now listed on the field's page, and each has a page of its own under its designation.

### Added
- **The field's bodies on its page and in the API (GEN.1105, part 2).** A field that keeps members lists them on its page under "The bodies it keeps" (designation, mass, diameter and shape, heaviest first), and each designation opens `/object/asteroid_field:<ID>/field-member/<rank>`: the body's shape, mass, diameter, density and makeup, its share of the field's mass, the history of the field it belongs to and a link to the heavier and the lighter body. The data is two new routes, `GET /api/objects/<ref>/field-members` and `GET /api/objects/<ref>/field-members/<rank>` (`describe/field_member_info.py`). Both are new routes, so the API version stays at 13. A free body has no star, so its page has no orbit.

generation-output: unchanged

TODO items resolved: GEN.1105 (with part 1: the members, their storage and designation, the pages and the API).
