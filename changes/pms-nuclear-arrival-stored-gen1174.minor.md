**Headline:** A star still contracting toward the main sequence now keeps both of its times, when its radius settles and when hydrogen burning arrives.

### Changed
- **The nuclear arrival age is stored (GEN.1174, schema v197).** `stars` and `bright_stars` gain `nuclear_arrival_age_gy`, the age at which MIST's tracks put a pre-main-sequence star on the main sequence. The star counts as pre-main-sequence only until its radius has settled (BHAC15), the earlier age `phase_end_age_gy` keeps, so the page of such a star can say both. The scatter writes the column, the fill that builds the system keeps it, the property sheet shows it ("Age it reaches hydrogen burning") only for a star that has one, and the uncharted-star API gives it. NULL for every star that is not pre-main-sequence and for rows made before v195; nothing is backfilled.

generation-output: changed
