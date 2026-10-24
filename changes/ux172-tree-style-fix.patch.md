**Headline:** The summary tree on the object pages is drawn properly again.

### Changed

- **Summary tree styling (UX.172).** A bad merge left a stray line and an unclosed block in the stylesheet, so the tree's rows, toggles and the rendered example only took their styles when the browser asked for reduced motion. The block is closed, and a test now fails for a conflict marker or an unclosed block in the stylesheet.
