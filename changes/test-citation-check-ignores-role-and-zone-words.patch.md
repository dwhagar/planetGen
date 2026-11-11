**Headline:** The citation check no longer reads "Manager, 2026" or "EST 2026" in the plan as an author and a year.

### Fixed
- **`NOT_NAMES` in `scripts/build_references.py` gains `Manager`, `EST`, `EDT` and `UTC`.** A role word or a time zone before a year ("Train Manager, 2026", "9:29 AM EST 2026") was read by the rewrite as a source named by author and year, and `test_no_citation_is_written_outside_the_registry` failed on main twice because of the TODO plan's wording. A new test in `test_citations_outside_the_registry.py` pins those four phrases. Reused: the existing `NOT_NAMES` stop list that the rewrite already consults; no new mechanism.

generation-output: unchanged
generation-affects: none
