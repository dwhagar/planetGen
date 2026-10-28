**Headline:** Every number an API request takes now states its unit in its name or its description, and a test keeps it that way.

### Changed
- **Units on API request fields (API.33, API part).** The request bodies' number fields that did not name a unit in their name (`position`, `num_orbits`, `type_id`, `sort_order`, `expires_days`, the worker `count`, the facility type fields' `min_value`, `max_value` and values) now carry a description starting "Unit:" in their schema, and the API reference has a section on the units of the numbers. `tests/test_api_field_units.py` fails over the schema of every body for a number with no unit.
