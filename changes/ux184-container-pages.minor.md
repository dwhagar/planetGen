**Headline:** Every globular cluster and superbubble now has a page of its own.

### Changed
- `/object/globular_cluster:<ID>` and `/object/superbubble:<ID>` show where a cluster or a superbubble lies in the galaxy, its stored properties (the radius of its sphere and where it was born), and an "About" text that says what it is, where its centre is, a cluster's mass, age and metallicity, and lists the star systems, nebulae and other objects inside it, each linked to its own page (UX.184).
- A globular cluster and a superbubble are now kinds an object reference can name: `globular_cluster:<ID>` and `superbubble:<ID>` resolve in `GET /api/objects/<ref>`, and their sheet and children (none) are served. A cluster's ID is never taken for a superbubble's.
- `GET /api/containers/<ID>` gives a cluster's mass, age and metallicity as `component`, and the contents list prints each object's own ID.
