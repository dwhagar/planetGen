### Changed
- **A lighter top bar with a settings gear.** The theme switch, search
  (when the bar has no room for it) and the account links moved under a
  gear at the upper right: Account, Admin, Generate, Stats and Logout
  for an admin, Login for a visitor (who never sees Stats). Galaxy,
  Sectors, Systems, Phenomena and Nav stay buttons while they fit and
  fold into a Menu when they don't, and the bar's search box shows only
  while its text entry is at least twice the width of its button. The
  layout follows the header's own width with container queries, so it
  needs no script; with script, an open menu closes on an outside click
  or Escape. The logged-in header no longer collapses at a wider width
  than a visitor's.
