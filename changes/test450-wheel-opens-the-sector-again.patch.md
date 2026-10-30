**Headline:** Zooming into a sector at the galaxy's core no longer gets stuck when the wheel points at a neighbouring empty cell.

### Changed
- **The wheel opens the sector again (TEST.450).** At the core the cells are small and close, so the wheel could land on a neighbouring cell nothing was generated in. The map asked the server for that cell's scene on every tick (a 404 each time), and each ask took the place of the real sector's scene load, which was then given up for good: the sector never opened and the zoom stopped about 1.4 parsecs out. An empty cell is now asked once, and a scene load that was only pushed aside is asked again.
