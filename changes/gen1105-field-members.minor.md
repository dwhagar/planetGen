**Headline:** An interstellar asteroid field made from a supernova's belt now keeps its heaviest bodies as individual members, each with a stable name made of the field's ID and its rank.

### Added
- **The stored members of a field (GEN.1105, part 1).** A field whose mass was solved with a ledger (the supernova host rule's) keeps the bodies the ledger stores, by the belt cut rule (every body over 1e19 kg and more of the next heaviest until 90 percent of the mass is held, at most 300), heaviest first, in the new table `field_members` (schema 222): rank, mass, radius, density, the four composition shares and the status word (`dwarf planet` where it is round, else `minor body`). A body takes the field's makeup and bulk density. `<field ID>.<rank hex>` names it (`object_uid.field_designation` and `parse_field_designation`, which refuse a belt's ID). A hand-made field and a Rogues top-up field keep none (the top-up field's largest body is under the cut). The field page, the member page and the API line (part 2) are not built yet.

generation-output: changed

TODO items resolved: none (GEN.1105 part 1 of 2: the members, their storage and their designation; the pages and the API version bump follow).
