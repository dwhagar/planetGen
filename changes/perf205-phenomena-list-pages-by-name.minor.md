**Headline:** The Phenomena list opens any page, sorted by name, as fast as the first, and every page of it is quicker to open.

### Changed
- **The Phenomena table pages by key when sorted by name** (its default). The page links carry the phenomenon the page before ended on, so a page deep in the list costs what the first does: the list walks each phenomenon table's name index and the best of those is merged in the database, where it sorted every table whole for each page. With 8,000,000 built phenomena the first page by name went from about 45 seconds to 10 milliseconds, and a page after a built row takes 4 to 6 milliseconds. `GET /api/phenomena` takes `after=<type>:<ID>` of a built row for it.
- **No page counts every phenomenon first.** A page that starts at the top, in any sort, no longer counts all the built phenomena to know where the scattered ones start; only a page that starts past some of them does. The other sorts (type, descriptor, radius, sector, on the Galaxy Map) still sort the whole list for a page.
- **Within equal names, a descending sort now lists the types and ids in reverse too.**
- **Two more name indexes** (asteroid fields and interstellar comets), built online by the update.
