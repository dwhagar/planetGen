**Headline:** The angles the region and window walks take for a small ball far from the core keep their digits, so a ball of a few parsecs at 262 kpc no longer loses cells to rounding.

### Fixed

- The walks that list the cells a ball reaches (`galaxy/region.py` through `spherical.sphere_shells` and `band_window`, and `window.cells_in_ball`) took each angle as the arccos of a cosine; below 1e-3 rad that cosine is within a few ulps of 1 and keeps hardly any digits (the angle of a 3 pc ball at 262 kpc was off by 4e-7 of itself). Below 1e-3 rad they now use the haversine form (GEN.376); larger angles are unchanged, and the cells a walk lists are still tested by distance, so no stored or generated value changes.

generation-output: unchanged
