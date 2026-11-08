**Headline:** The Chart View and side strip browser tests open the Chart View sub-menu and expect the isometric slant.

### Fixed
- **Chart View browser tests (TEST.696).** The Chart View and side strip controls now sit in a sub-menu that starts closed, so twelve browser tests waited two minutes for a button that was not visible; they now open it first. The plane view test expected a 62 degree slant where the Chart View is isometric (54.74 degrees), and one test read the plane's height before the grid had been rebuilt; both are corrected.

generation-output: unchanged

TODO items resolved: TEST.696 (plane view and side strip).
