**Headline:** The search, nearby and chart forms and the table filters use one field look, with the server's error shown on the field it is about.

### Changed
- **One form-field look and error slot (UX.124).** A new `--line` colour (3:1 or more against the panel in light and dark) is the edge of every field, and the macros in `partials/ui_fields.html` draw a labelled Shoelace `sl-input`, `sl-select` or `sl-checkbox` with the server's validation error under it. The search page, What's nearby, the Plot Course chart form and the data tables' "Starting at" box and filter menus use them; the What's nearby distance error now shows under the distance box. With script off each field is a plain control in a `<noscript>`, so every form still posts. Fields with a `datalist`, and the hidden, submit, password and one-time-code fields, stay native.

TODO items resolved: none. Added: none.
