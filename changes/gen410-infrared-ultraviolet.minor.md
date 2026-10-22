**Headline:** The gas table now holds each gas's infrared bands and its ultraviolet and visible absorption.

### Changed
- **Infrared bands and ultraviolet bins in the gas table (GEN.410, third slice).** 37 gases carry their infrared band centres (part 03 of the atmosphere research) with the saturated-absorption rank the research states, and 22 carry their mean photoabsorption cross section over 35 bins from 10 to 1100 nm, with where the data stop, where they come from and the limits the research names (the structured O2 bands, the visible methane bins). `gas("O3").uv.at(255.0)` gives the Hartley band. Nothing reads them yet.

generation-output: unchanged
