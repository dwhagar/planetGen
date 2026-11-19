**Headline:** The far-from-center buffer test holds the whole-galaxy block prisms to the float precision rule that the 10 kpc limit stands for.

### Fixed

- `test_layers_far_from_the_center_are_built_relative_to_origins_of_their_own` held every layer to `SAFE_EXTENT_PC` (10 kpc), including the "block prisms" of the one whole-galaxy view, whose two meshes (solid and glass, MAP.150) are 15,540 pc across. The limit is a float32 precision bound, and `camerarelative.js` says that view draws its blocks as one part that spans more, where a pixel is 100 pc or more. The test now checks that layer by the rule itself: the float32 step at the part's farthest vertex (0.001 pc here) stays under a tenth of the 100 pc pixel. Every other layer keeps the 10 kpc check, and an opened sector's blocks are still held to it. The two-mesh design is kept (MAP.432).

generation-output: unchanged
