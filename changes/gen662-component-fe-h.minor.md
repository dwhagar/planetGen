**Headline:** The stars of a globular cluster, a young cluster or a satellite galaxy now draw their composition from their own component instead of a halo value.

### Changed

- **Cluster and satellite stars take their component's [Fe/H] (GEN.662).** The stars the component pass places in a cluster or a satellite draw their birth [Fe/H] from that component's own value (with the scatter of its family) before their lifetimes and giant branch are worked out, and the row stores it. A star with no component given takes the generic cluster or satellite distribution; the halo's value, and the reference 0 their trees were read at, are gone from those paths.

generation-output: changed
