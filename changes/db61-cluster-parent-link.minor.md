**Headline:** A cluster can now be unnamed, and a star can point to the cluster that is its parent, ready for the deep-sky work.

### Changed
- **Unnamed clusters and `cluster_id` on `stars` and `bright_stars` (DB.61, schema v127).** `galaxy_components.name` may be empty, with the unique key kept over the named rows, and each star table gains an indexed `cluster_id` that refers to a cluster. Existing clusters keep their names and every star's link starts empty. The migration rebuilds the two star tables once (the time grows with their rows); nothing writes or reads the link yet.

TODO items resolved / Added: DB.61 / Added: none
