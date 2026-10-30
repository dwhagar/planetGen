**Headline:** The data behind an asteroid belt information page now exists: the belt's individually tracked bodies by designation, what the rest of its mass is made of by size, and how thick the belt is in objects per square astronomical unit.

### Added
- **`physics/belt_page.py` (UX.240).** One plain structure with three blocks over a belt's ledger: the stored bodies (designation, mass, diameter, distance, taxon and composition, density, status word and Soter's mu), the rest of the belt mass as diameter bins (count, mass and share, whether drawn on demand or a total only, typical composition and density) that add up with the listed bodies to the belt's mass, and the density factor (objects per square AU above 1, 10 and 100 km, and kilograms per square AU). Every number carries its unit. Nothing serves it yet; the page route and the API field are UX.241.
- **`belt_member.typical_density_kg_m3`**: the bulk density of a typical body of a taxon at a mass, for the size table.
- **Stand-in rows BP1 to BP4** in the gap-filling plan (section 6.120), the stand-in list and the wiring register.

generation-output: unchanged
