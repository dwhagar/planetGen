**Headline:** The map menu test accepts the View menu's group names and the Gravity note as text.

### Fixed
- **System menu browser test (TEST.696).** The check that every visible word in a map menu is a heading or a control knew only the old heading class. The System Map's View menu now names its groups with `view-menu-group-name` and carries a one-line Gravity note, so the test flagged them. It accepts both and says the element's class and parent when it flags one.

generation-output: unchanged

TODO items resolved: TEST.696 (system menu).
