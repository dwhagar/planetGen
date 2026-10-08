### Changed

- Choosing a NAV course on the Galaxy Map and the Sector Map stays on the map
  (NAV.29, NAV.33). An object's panel offers "Start Here" and "End Here"
  instead of "Nav from here" and "Nav to here". The first one pressed keeps
  the view and zoom where they are, puts the chosen end in a banner and in the
  URL (`?pick=to&from=system:12`), and asks for the other end; the user zooms
  out, in and across to find it with the same picker, and its button opens the
  NAV page with both ends. Cancel on the banner clears a pick begun on the map.
  The same buttons end a pick the NAV page began, and the system and
  phenomenon pages' pick buttons carry the same words. The sector scene JSON
  no longer carries NAV links: the page's script builds them (`navpick.js`).
