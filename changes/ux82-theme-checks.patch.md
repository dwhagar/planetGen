**Headline:** The Shoelace primary-button hover and the dark theme's subtle background now meet the theme checks carried over from PR #800.

### Fixed
- **UX.82: the primary button's hover fill was 4.23:1 against its label in the light theme, and `--bg-subtle` was a colour with alpha in the dark theme.** The hover fill (`--sl-color-primary-500`) now mixes the accent with the text colour instead of the panel, which clears 4.5:1 in both themes; the placeholder (4.41 before #800) already measures above 6:1 in both. `--bg-subtle` in the dark blocks is the opaque blend of the old colour over the page, so what it sits on no longer changes its shade. The map viewports' token block also gained the two references tokens it lacked, which fixes `test_map_black.py`, red on main since the references colours were added.
- **New checks (`test_theme_checks.py`):** every icon file draws in `currentColor`, both contrast figures hold, and `--bg-subtle` has no alpha. Nothing generated changes.
