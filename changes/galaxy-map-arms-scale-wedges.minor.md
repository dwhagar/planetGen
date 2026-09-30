### Added

- Galaxy Map: wedge lines run out from the galactic core in the plane, each
  labelled with its bearing (degrees counterclockwise from +X, ring slot 0),
  with a Wedges button to hide them.
- Galaxy Map: the scale readout has three lines, what one screen pixel
  spans, how big one block is, and a bar, each in sectors, pc and ly.

### Changed

- Galaxy Map: the density prisms are shaded by each prism's arm factor
  (its density over the ring's mean) as well as its density, so the spiral
  arms stand out at every zoom.
- Galaxy Map: the density blocks fill their whole cells, so the galaxy is one
  solid made of blocks with no gaps. A Slice button (on by default) cuts the
  solid at the focus's layer, so the view looks down on its cut face;
  turning it off shows the whole solid.

### Removed

- The server's leftover density point clouds: `/api/galaxy/tiles` and
  `/galaxy/tiles` no longer take `density=` or return `density`, and
  `galaxyViewport.density_sample_points` / `density_points_for_tile` are
  gone. The map has drawn density itself since the prisms.
