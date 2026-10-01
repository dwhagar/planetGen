### Fixed
- **Buttons always have space between them (UX.16).** One shared rule in
  `static/style.css` gives every group of buttons on every page (`.btn`,
  `.btn-small`, `.starmap-btn`, plain `<button>`s, and side-by-side forms
  that each hold one) the same gap, across and between wrapped lines:
  `--btn-gap`, 0.5rem with a mouse and 0.75rem on touch screens
  (`pointer: coarse`), so 44-48 px touch targets never sit edge to edge.
  Groups that had their own smaller gap (the Galaxy Map's address
  matches, 0.3rem) now use it too.
- **An object's data sits beside its 3D render when there's room
  (UX.15).** On the System Map the info panel moves to the right of the
  map once the panel is at least 46rem wide (the map shrinks to fit and
  stays square); on a phenomenon page with a 3D view (neutron stars,
  black holes, quasars, rogue planets, interstellar comets) the data
  table sits beside the view from 50rem. Narrower screens keep the
  stacked layout. Both switches are container queries, so the layout is
  set before the render loads and nothing jumps.
