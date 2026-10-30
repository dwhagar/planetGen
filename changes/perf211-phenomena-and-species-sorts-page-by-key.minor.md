**Headline:** Every sort of the Phenomena table and the Species table now answers a page in milliseconds, however many rows there are.

### Changed
- **The Phenomena table** sorted by type, descriptor, radius, placed or sector used to sort every phenomenon of all eight kinds for each page (the first page by type took 56 s with 8,000,000 phenomena). Each kind now has indexes for these sorts, built when you run ./update.sh (migration 0186: expression indexes make the database hold writes to that table while they build, so run it with generation stopped), and the first page and every page after a row take 3 to 41 ms. The sector sort reads the sectors by name and sorts only their phenomena. Where two rows tie in a sort, the order now follows the direction of the sort (name, then type, then ID); before, a descending sort listed ties from A to Z.
- **The Species table** sorted by homeworld or system now walks the planets or systems by name and sorts only their species (1,000,000 species: 8 to 13 ms for a page deep in the list).
- **`GET /api/phenomena?after=<type>:<ID>`** carries on in the order of `sort`, not only by name; the table's Next link does so in every sort.
- **`planetgen sectors` and `planetgen systems`** list 100 rows at a time like `planets` and `moons`: `--limit N` (0 lists all) and `--after ID`, and the last line prints the `--after` to carry on with.
