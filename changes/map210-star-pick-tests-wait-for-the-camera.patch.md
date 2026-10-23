**Headline:** The star selection browser tests wait for the camera to stop and also cover a sector 8 kpc from the core.

### Changed

- The opened-sector tests now wait until the camera holds still before clicking (a click on a moving map is ignored, and a flight to a far sector lasts longer), work out each star's screen position again at the click, and click only stars inside the window.
- New test: a star of a sector 8,000 pc from the core is selected by a click, the panel naming it and linking its system.
- Boss's report that nothing could be selected in sector view did not reproduce on current main: every star of the opened sector picked correctly in twelve fresh test databases, near the core and far from it.

TODO items resolved / Added: TEST.281 / Added: none
