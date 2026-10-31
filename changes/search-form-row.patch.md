**Headline:** The Search page's Name box, Stars selector and Search button now sit together on one row.

### Fixed
- **The Search page left a large empty gap above Stars and Search.** The Name box's field wrapper carried a row-sized `flex-basis` that acted as a 12rem minimum height inside its column, stretching the row and pushing Stars and Search to the bottom right. The wrapper now sizes to its content, and Stars keeps a fixed width beside it; the row still wraps on narrow screens.
