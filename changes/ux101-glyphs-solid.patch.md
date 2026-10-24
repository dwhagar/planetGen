**Headline:** The habitability symbols are solid colour-shaded shapes you can read in both themes, and the radiation symbol is a trefoil.

### Fixed
- **Every icon of the sprite drew as a solid black shape.** The sprite set its outline in a style element of its own, which the site's Content-Security-Policy (`default-src 'self'`) refuses, so the shapes took the default black fill and the PHI-6 symbols showed only a tiny coloured dot (Boss, 9:07 PM EDT 2026-10-23). The outline is now the `.icon` rule of `style.css`.
- **The PHI-6 symbols are solid shapes with no border (Boss, 9:10 PM EDT).** Each is filled in its rating colour with its inner details cut out, at 0.8 opacity so the page behind shades the colour in the light and the dark theme, and a little larger (1.45em). The accessible labels are unchanged.
- **The four life symbols are solid and borderless too (UX.151).** Microbial, Flora, Fauna and Sentient are filled shapes at 0.8 opacity like the PHI-6 ones; the number beside each is unchanged.
- **The radiation symbol is a trefoil.** A dot and three blades replace the three arcs that could not be read as radiation.
