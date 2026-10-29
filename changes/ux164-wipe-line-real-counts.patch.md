**Headline:** The line after a reset now says how many tables were wiped and how many rows they held, counted, instead of "about 1 row".

### Fixed
- **UX.164: the wipe line gave the database's stale statistic ("Wiped 75 tables (about 1 row)").** `planetgen reset` and the Generate page's reset now count each table before wiping it (a table the statistics put under 200,000 rows is counted exactly in a moment) and say, for example, "Wiped 75 tables from 'galaxy' (12 tables held rows, 4,512 rows removed)", or "all were already empty" for an empty database. The few very large tables keep the statistic, as counting them row by row would take minutes, and the line says so ("... and about 9,000,000 more by the database's statistics in its large tables"). No stored value changes.
