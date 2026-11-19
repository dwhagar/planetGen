**Headline:** A New galaxy job no longer stops at its luminosity scatter, which took the backfill's sky stars for a luminosity scatter already stored.

### Fixed
- **The luminosity scatter of a job refused the job's own backfill (OPS.280).** In a job where a backfill comes before the massive-star scatter, the massive-star scatter already refuses a galaxy that holds a scatter and empties the table before the backfill; the luminosity scatter after it then found the backfill's sky stars under the mass limit and stopped with "The bright-star scatter is already stored". A luminosity scatter prepared that way no longer checks. A luminosity scatter on its own still refuses a stored one. A test runs the job's chain end to end.
