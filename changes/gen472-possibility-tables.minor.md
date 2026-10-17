**Headline:** Every composition table now loads through one gated possibility-table format, and a draw takes the rows that pass the body's class and conditions and picks one by normalised weight.

### Added
- **The gated possibility table and its loader (GEN.472).** A row names a possibility, the classes it supports, optional conditions on stored quantities, a weight (a number, `open`, or a function of the star and the formation distance), its source with the verification state, and its licence. A draw drops the rows that fail the gate or weigh zero, normalises the rest to 100 percent and picks one on the generator's own random stream. A class no row supports, or a table whose weights are all zero, is an error naming the class and the table, and a row without a source or a licence is refused when the table loads. Open weights share what the sourced weights leave equally and are listed with the stand-in rows in one report. Fractions inside a row are taken directly from a sourced kernel or drawn from a Dirichlet distribution where only a range is known, and every component must be a row of the one shared materials table.
- **`draw.gamma` and `draw.dirichlet`**, built on the generator's own random wrapper like the other draws.

TODO items resolved: GEN.472. Added: none.
