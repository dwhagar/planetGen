**Headline:** A generated galaxy now has young open clusters, drawn from its own star formation, with stars of their own age.

### Added
- **Young open clusters (GEN.244).** The disk's star formation history (GEN.469) says how much mass formed in each ring in the last gigayear; the cluster formation efficiency of the ring's star formation rate density (Kruijssen 2012) puts a share of it in bound clusters, with initial masses falling as M^-2 from 50 to a million solar masses (Lamers et al. 2005). Stellar evolution and the tide of the place take mass from each (Lamers, Gieles and Portegies Zwart 2005), and one that is gone is not drawn. A cluster lies in the disk at the radius of the ring that formed it, on the arms when young, on the circular orbit, with the disk's [Fe/H] there, and a half-mass radius that grows five to ten times over its first 100 Myr (Portegies Zwart et al. 2010). Nothing is taken from the Milky Way's catalogue: a drawn galaxy has from a few thousand to tens of thousands, and the most compact have the most.
- **Twelve age-bin populations.** All the stars of a cluster share its age, so they belong to the population of the quarter-decade of age it falls in (`open_cluster_00`, 1 to 1.78 Myr, to `open_cluster_11`, 562 Myr to 1 Gyr). The cluster adds its density to its cell in that population, the profile and the sector page name it ("Open cluster stars, 1 Myr to 1.78 Myr"), and the component pass draws its bright and massive stars. OB associations stay the disk's young population.

### Changed
- `bright_stars.population` takes the twelve new names (migration 0112, no row is touched).
- The plan draws and stores the young open clusters with the globular clusters, so a new galaxy's stars, systems and sector profiles differ near them.
- The galaxy's component field is built again only when the stored components' signature changes (their count, first and last id and mass sum), so a page does not read tens of thousands of cluster rows on every request.

generation-output: changed
