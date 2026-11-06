**Headline:** The System Map's and the sector page's Menu now fold into headings with counts, real checkboxes and a mark where a switch is not as it starts.

### Added
- **Shared View menu builder (UX.278).** `static/viewmenu.js` makes the Menu of a map from the layers it has, using the headings each layer
  declares (UX.277): a fold per heading (Show, Color and fill, Chart and side views, Move and Settings; a map shows only the ones it has
  layers for), a count in each ("Show 4 of 12 on"), real checkboxes that keep the layers' `data-layer` and `data-*` attributes, an
  all/some/none box for a group of two or more layers, the word "changed" beside a switch that is not as it starts, and the open or folded
  state of each heading kept in the browser. A menu of twelve layers or fewer starts open; a longer one folds all but the first heading.
  The System Map and the sector page use it; the Galaxy Map's Menu is unchanged until UX.280. The search over the entries is UX.279.
