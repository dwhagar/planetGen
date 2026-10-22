**Headline:** Bookmarks now sit on a shared saved-list store, ready for saved courses and for lists kept in an account.

### Changed

- The Bookmarks code is split in two: `savedlist.js`, a store for one list of saved entries behind a backend (the browser's own storage today, an account later), and `bookmarks.js`, the star button, the Bookmarks menu and the NAV page's select on top of it. Saving, renaming and deleting a bookmark look and behave as before, and your existing bookmarks are kept (the same browser storage key and format).
- A second kind of list, such as saved courses, needs no new storage code.

TODO items resolved / Added: UX.129 / Added: none
