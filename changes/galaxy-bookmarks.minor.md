### Added

- Bookmarks (MAP.23, which also finishes the NAV page's map picks,
  MAP.22, and so the Galaxy Map drill-down, MAP.2). A ☆ on the Galaxy
  Map's breadcrumb saves the view (its stage URL) or the selected sector,
  and a ☆ Bookmark button on system, phenomenon and sector pages saves
  that page; it turns into ★ Bookmarked, and pressing it again offers to
  remove the bookmark. A Bookmarks menu in the map's controls opens,
  renames and deletes them, and Ctrl+1 to Ctrl+9 open the first nine
  while the map page has the focus (not while typing in a box). The NAV
  page offers the system and phenomenon bookmarks as a start or a
  destination, and a sector bookmark opens that sector's system picker.
  Bookmarks are kept in this browser (`localStorage`, one list per
  database, up to 100), in the new `static/bookmarks.js`; with storage
  blocked there are just none.
