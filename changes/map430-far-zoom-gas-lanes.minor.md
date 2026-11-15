**Headline:** The Galaxy Map at galaxy scale now draws smooth gas and dust lanes from the galaxy model's gas density (the arms and the inner disk), which fade as the camera comes in and the nebulae of the nodes take over.

### Added

- MAP.430: one more volume on the Galaxy Map, the galaxy's gas on the plane (`GET /galaxy/nebula-volume?key=gas`, `db/nebula_volumes.gas_lanes_volume`): a 128 pc lattice sampled from the same model gas every nebula contour reads (`nebula_field.gas_factor_array`, GEN.1267), bent and clumped by a smooth warp and noise, drawn in the dark family's colour through the nebula volume layer's capped sum. Each lane splat fades by its width on the screen, so the lanes give way to the nodes' nebulae as the camera comes in and are never a grid of big blobs. Display only.

generation-output: unchanged
generation-affects: none
