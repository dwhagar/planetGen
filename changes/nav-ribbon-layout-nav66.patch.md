**Headline:** The schematic route map has its layout: a route and the links round it become a ribbon of stops, hops and stubs.

### Changed

- New `planetgen.web.maps.navribbon.layout_ribbon(route, width)` (NAV.66, part 2 of the schematic route map): pure arithmetic over a route as `nav_between` and `GET /api/nav` give it, with its `graph` when it has one. It returns the stops left to right with their kind (a system, a stop in unknown space by its kind, or a phenomenon end) and role, the hops with their length written on them (the site's unit text and its first unit alone) and the flag for unknown space, the neighbour links from `route.graph` as stubs above and below each stop (alternating, centred, with an "and N more" text where the graph's cap cut some), and the links between two stops that are no hop as arcs.
- A hop's gap on the page follows the logarithm of its length (68 px at 1 ly or less, 200 px at 10,000 ly or more), the same on every route, so one 6,500 ly hop no longer squeezes nine hops of 3 ly to dots. A ribbon wraps to a new row where the next stop would pass the right margin of the panel (at least 240 px wide); a hop that wraps is two pieces, one to the right edge and one in from the left edge, and says so. The same route gives the same answer, with coordinates rounded to a hundredth.
- Nothing draws it yet: NAV.67 puts it on the NAV page as an SVG. Tests on fixed routes: the gap law, a plain route, stops in unknown space, no route and one stop, wrapping and narrow panels, stubs, the cap's text, arcs and links across rows, determinism, and a route from `nav_between`.

generation-output: unchanged
