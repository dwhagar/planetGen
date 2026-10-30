**Headline:** `planetgen planets` and `planetgen moons` show a page at a time instead of reading every planet or moon in the galaxy.

### Changed
- **The planets and moons listings take `--limit` and `--after`.** Each prints 100 rows by default and, when more match, the `--after` to carry on with; the next page starts after that row in the same order, so a page deep in the list costs what the first does. Before, the commands read every matching row, which at a full galaxy meant millions. Filtered by system, sector or class a page takes about 40 milliseconds with 6,000,000 planets. Unfiltered, the order (system name, then orbit) still sorts every planet, so the first page of the whole galaxy takes seconds.
