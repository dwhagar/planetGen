**Headline:** The text description of a body now keeps each quantity in the unit it is stored in, and converts it only when it is shown.

### Changed
- **Describe quantities** (GEN.574): a quantity in the description model is `{"kind", "value", "unit"}` instead of `{"kind", "si"}`. `model.qty` checks the unit against the unit registry and `model.measured` keeps the stored quantity's own unit, so a length held in megametres is no longer converted to metres on the way in. The wiki, text and web renderers convert to the kind's base unit when they show it; the page's `data-si` attribute, which the browser reads, is unchanged. What you see on a body page does not change.
