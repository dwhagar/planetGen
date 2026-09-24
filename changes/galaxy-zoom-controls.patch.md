### Fixed
- **Galaxy Map clicks recentered the view far from where you clicked.**
  Clicking empty space picked a point on a sphere the camera itself sits
  on, so the view jumped thousands of parsecs away (often out of the
  galaxy entirely). It now lands on the spot under the cursor in the
  galactic plane (or facing the camera when the disk is seen edge-on),
  kept inside the galaxy.
- **Double-click zoomed somewhere other than where you double-clicked.**
  Each of its two clicks recentered again before the zoom; now only the
  first click recenters and the double-click zooms in on that point.
- **Clicking a sector centered on the edge of its marker, not the sector,**
  so it slid off-center and out of view as you zoomed in. It now centers
  on the sector's own position.
- **Scroll zoom ignored how far the wheel moved.** A trackpad's stream of
  tiny scroll events each took a full step, making zoom race; zoom now
  follows the scroll amount (about 1.28x per mouse-wheel notch, pinch
  supported).
- **Zooming all the way out didn't show the whole galaxy.** The farthest
  zoom now backs the camera off far enough for the galaxy's whole disk to
  fit in the map.
