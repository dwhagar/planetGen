**Headline:** The one-off system page can now generate a rogue planet or an interstellar comet, not only a star system.

### Changed
- **A Result choice on the one-off system page (UX.298).** Under Options, Result picks a star system (the default), a rogue planet or an interstellar comet. The two starless results come from the existing `planetgen phenomenon` generator, which now takes the same `--output` and `--format` that `planetgen system` has: it writes the description and saves nothing. Only the name and the debug log apply to a starless result, and any other filled-in option is refused with its name instead of being ignored. The page shows it as it shows a system, with the same downloads.

generation-output: unchanged
