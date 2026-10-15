**Headline:** Every multi-star and star-like-object system type in the catalogue now has an entry in one rate table, ready for the member-tree generator to read.

### Changed
- A new table, `SYSTEM_TYPE_RATES`, holds one entry for each of the 57 types in the system-type catalogue (stellar architecture, interacting binaries, white-dwarf systems, compact-object systems, substellar companions, runaways): the reference class, default, low and high, evidence tag, hosting populations, generator, scatter or fill mode, theoretical flag and audit numbers. The tuning knobs for the tight-period tail, post-common-envelope share, accreting share and merger-product shares sit beside it. Nothing the generator makes changes yet; a default is now an edit to this table.
- A test reads the catalogue and fails when a row has no entry, an entry lacks a generator, owner or audit, a range is out of order, or a number with no value of its own does not say where it comes from.

TODO items resolved: GEN.268 / Added: none
