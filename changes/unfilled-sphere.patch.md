### Fixed
- **Zoomed in on the Galaxy Map, unfilled sectors still drew as a ball.**
  Two things drew it. Unfilled (not-yet-generated) sector dots were
  clipped to a 20 pc ball around the view center while the view reached
  out to 200 pc, and near the galactic plane every slot qualifies, so the
  ball was solid. They now fill the whole view, shown while the view
  radius is 32 pc or less, and fade out toward its edge. Wider zoomed-in
  views get the density cloud instead, which had its own ball: a view a
  few hundred parsecs across kept almost none of its galaxy-wide samples
  and topped up with uniform points around the view. It now samples such
  views locally against the real density, so the cloud follows the disk.
