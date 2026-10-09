### Changed
- System and wiki Markdown is rendered by the `markdown` package (cut down to headers, pipe tables, paragraphs and `<sup>`), replacing the hand-written `mdconvert.py`. The pages look the same; table alignment colons now give an `align` attribute instead of being ignored (UX.39).
