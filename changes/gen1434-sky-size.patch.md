**Headline:** the sheet of a planet or a moon now shows how wide its sun looks in its sky, a moon's sheet how wide its planet looks, and the sheet of a body in a multiple system how wide each other sun looks.

### Changed

- `describe/sky_size_rows.py` (GEN.1434, new) and `describe/sheet.py`: derived rows built when the sheet is made, no stored column and no migration: "Width of its sun in the sky" (planet and moon; the angular diameter `2 asin(R / d)` at the circle that gives the same mean light as the orbit, 0.533 degrees for the Sun from the Earth), "Width of its planet in the sky" (moon; the Earth is 1.9 degrees wide from the Moon, 3.7 times the Moon's width from the Earth) and up to four "Width of its second sun in the sky" rows for a body whose system has other sources bright enough to show (the same list and widths as the sky text of GEN.1443, `sky_view.bodies`, read from the loaded system only for a body that has stored other suns). A rogue planet has no sun and no row.
- Registry rows SZ1 and SZ2 added (section 6.280 of the gap plan, provisional until the train numbers it).

generation-affects: none
