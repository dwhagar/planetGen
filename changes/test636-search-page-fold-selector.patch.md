**Headline:** The search page form test opens the fold it means, not both.

### Fixed
- **Search page fields test (TEST.636).** `test_each_kind_of_field_is_filled_and_submitted` clicked `details.search-more > summary`, which has matched two folds ("Search by object and size" and "Search by physical values") since the physical-value boxes came in. It now opens the size fold, which holds the boxes the test fills. No change to the page.
