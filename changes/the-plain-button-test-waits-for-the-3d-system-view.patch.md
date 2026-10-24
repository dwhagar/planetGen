**Headline:** The every-button browser test waits for the 3D system view to finish loading before it clicks.

### Fixed
- **`test_every_plain_button_does_something[web.system]` timed out on a body-list button of the 3D view (TEST.350).** The system page opens in 3D by default and rebuilds its "Bodies in this view" list once the scene has loaded, so a button clicked just before that was "not stable" and the click waited out its two minutes. The test now waits for the 3D view to be ready (or hidden) after every page load.
TODO items resolved: TEST.350. Added: none.
