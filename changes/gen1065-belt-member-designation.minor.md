**Headline:** A member of an asteroid belt now has a designation, the belt's ID, a point and the member's rank in hex, that is a pure function with a round trip.

### Added
- **Belt member designation (GEN.1065).** `galaxy/object_uid.designation(belt_id, rank)` prints `<belt ID>.<rank hex>` (rank 1 is the heaviest member by the mass drawn at birth; at least four hex digits), `parse_designation(text)` reads it back through `parse_id`, and `is_on_demand(rank, individual_count)` says whether a rank names a body drawn on demand and never stored. It reads no stream and no naming key, so a member has no name, only its designation, and one text. Nothing is stored or shown yet; the display and storage are the lock items of the belt evolution chain.

No generator output changes.

generation-output: unchanged

TODO items resolved: GEN.1065.
