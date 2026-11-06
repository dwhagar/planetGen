**Headline:** The stars of a cluster now point to their cluster, and a cluster gets its name when the first sector that holds one of its stars is filled.

### Changed
- When the plan places the stars of a cluster, each star now carries the cluster it belongs to, and the system built around that star when its sector is filled passes the cluster to all of its stars. Field stars and the stars of satellite galaxies carry none, so a cluster page and the maps can list the stars that belong to it.
- A generated cluster has no name until then. The first fill that touches one gives it a name, once: two words from its Object ID under the galaxy's naming key, spelled with its sector if another cluster already has that name. A second worker or a resumed fill finds the name already there and changes nothing. A cluster with a catalogue name keeps it. Until a cluster is named its page and labels show "Globular cluster" and its ID. The stars drawn for each cluster are drawn from a stream of the cluster's own seed instead of the old placeholder name, so a galaxy planned before this change draws its clusters' stars differently.

generation-output: changed
