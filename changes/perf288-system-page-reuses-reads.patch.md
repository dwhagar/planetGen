**Headline:** A system page sends two fewer database queries, because it reuses the star rows and the belt rows it has already read.

### Fixed
- `system_detail` read the primary star's heliosphere radius with a query of its own although its star rows (read a few lines before) already carry it, and read the asteroid belts a second time to number the irregular moons' designations. Both reuse the rows already read now, so a system page sends at most 19 statements with one host (it was 20 or 21 on main, over the 19 that `test_sector_and_system_pages_use_a_fixed_number_of_queries` allows). The page shows the same data, and the limit is unchanged.
