**Headline:** A whole-sphere galaxy scatter now runs to the last shell instead of stopping at 97 percent.

### Fixed
- **The massive-star scatter stopped at 97 percent with "Out of range value for column 'band_index'" (DB.72).** A shell's bands run from `-(h + 1)` to `h + 1`, which passes the 32,767 that a SMALLINT holds at shell 20,861 and reaches 102,942 at shell 65,535, the outer shell of the 262 kpc reach. `band_index` is now INT in `bright_stars`, `sector_stats` and `phenomenon_scatter` (migration 0108, schema v108; each table is rebuilt in place and no value changes). A test stores and reads back a row at the outermost shell's polar cap band in each table, and another refuses any column named for a band that is narrower than INT.

TODO items resolved: DB.72. Added: none.
