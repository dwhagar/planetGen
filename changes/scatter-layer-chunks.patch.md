### Changed
- The plan's scatter passes now hand their layers to the workers in a few chunks instead of one task per layer, so two or more workers no longer spend more time starting tasks than drawing stars. Rows are unchanged.
