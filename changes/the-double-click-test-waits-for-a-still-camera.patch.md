**Headline:** The double-click flight tests wait for the camera to hold still, and the slow-first-click test says what happened when it fails.

### Fixed
- **`test_a_double_click_flies_even_when_the_page_is_slow_to_answer_the_first_click` failed once under load (TEST.262).** It passed alone and in 41 runs under a four-way browser load, so the cause was not found. The tests of the flight to a star now wait until no flight is under way and the camera has not moved for a quarter second before they aim at the star, which `_ready` did not guarantee, and the slow-click test's failure message now carries the times of the pointer events and the camera before the click.

TODO items resolved: none. Added: none.
