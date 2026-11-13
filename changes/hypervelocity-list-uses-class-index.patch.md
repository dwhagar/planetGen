**Headline:** The Lists page and the cache warm-up no longer time out on the hypervelocity star list, because a scattered class is found through its own index instead of a walk of the whole table.

### Fixed
- **`/lists` answered 504 and stopped the cache warm-up.** The list of hypervelocity stars (and every other uncharted kind) picked its rows from `phenomenon_scatter` with `COALESCE(subtype, '') = ?` over the OR of the classes. A function of a column cannot use `idx_phenomenon_scatter_class`, so MySQL read the table in id order from the start, and a class whose rows sit late in it ran past the 10 second statement limit. The page of a kind now asks one class at a time on the columns themselves (`kind = ?` and `subtype = ?` or `IS NULL`), and the Phenomena table, its count and the search panel use the same clause (`query.scatter_class_clause`). On a 3,000,000-row copy the hypervelocity page went from 1.3 s with a full walk to 0.03 s on the class index.

generation-affects: none
