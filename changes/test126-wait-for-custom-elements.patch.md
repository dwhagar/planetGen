**Headline:** The "every plain button does something" browser test waits for the page's custom elements, so it no longer clicks a button that vanishes as its dialog is defined.

### Changed

- `test_every_plain_button_does_something` took the page's buttons as soon as the page loaded. Until the page's custom elements are defined, the children of a shut `sl-dialog` are drawn in the page, so a dialog's "Close" button counted as visible and then disappeared when the dialog was defined, leaving the click waiting out its 120 second timeout (seen on the sector page, 2 of 6 runs alone; TEST.126). The test's settle step now also waits until no custom element is left undefined. Three runs alone pass.
