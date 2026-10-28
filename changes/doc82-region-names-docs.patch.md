**Headline:** The named regions now have their own design page and a check that neighbouring regions never share a name.

### Changed
- **Region names documented and audited (DOC.82).** A new page, `docs/design/galaxy/named-regions.md`, says how a region is named, framed and stored (the measured numbers behind the choices included), and the pages for the block ladder, object pages and the map wire format point to it. The release audits gain a check that no two neighbouring regions share a full name, a first word across two families or an ordinal inside one (`tests/audit/region_names.py`, with its own tests). Nothing visible changes and no generated value changes.
