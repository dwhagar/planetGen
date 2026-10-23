**Headline:** The number of globular clusters a galaxy holds now follows its halo mass instead of always being the Milky Way's 157.

### Changed

- Every galaxy drew exactly 157 globular clusters, whatever its size. The count is now the cluster system's mass over the mean cluster mass: the system is a fraction `eta` of the halo mass, drawn per galaxy log-uniform between 2.9e-5 and 6e-5 (Harris 2017, Hudson 2014 and Harris 2013 disagree by a factor of two) with a 0.2 dex scatter (Hudson 2014), and the mean cluster mass comes from the generator's own truncated log-normal mass function. A Milky Way mass halo gives about 150 clusters, a tenth of that mass gives about 15, and a halo too small for one cluster gives one or none by a random rounding that keeps the mean exact. A galaxy without a halo has no clusters. The count law takes a family factor (1 for now) so a later galaxy family can scale it. `tuning.CLUSTER_COUNT_MILKY_WAY` is removed. The clusters' own masses, sizes and places are unchanged.

generation-output: changed

TODO items resolved / Added: GEN.403 / Added: none
