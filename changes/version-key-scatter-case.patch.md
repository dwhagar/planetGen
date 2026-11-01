**Headline:** The version-key step of an update no longer spends minutes building bounds nobody reads.

### Fixed
- **The bright-star scatter test case of the version key took minutes (case 5 of 8 in the update's progress bar).**
  The share of a population's stars that are bright was read from a full sampling table, whose per-cell bounds are most
  of the cost of a table and are never used for a share. It is now worked out from the quadrature alone
  (`star_population._bright_fraction`), to the same numbers, so every case digest and every version key is unchanged.
