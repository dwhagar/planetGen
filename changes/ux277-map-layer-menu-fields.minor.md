**Headline:** Every map layer now declares the View menu heading and search words it belongs under, ready for the grouped menu.

### Added
- **Layer registry fields for the View menu (UX.277).** A layer of the Galaxy, Sector and System Maps (`static/maplayers.js`) now names its
  `heading` (Show, Color and fill, Chart and side views, or Move and Settings), its `keywords`, a `common` flag for the few most used
  switches and an optional `settings` hook for a layer that needs a slider or a filter. Every existing layer has its heading and words,
  and `registry.unheaded()` names any layer without them, which browser tests on the three maps hold empty, so a layer added later cannot
  skip them. `registry.entries()` hands the menu builder (UX.278) and the search (UX.279) what they need. No look change.
