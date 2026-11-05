**Headline:** The NAV page's course diagram now shows each stop's nearest neighbours as faint stubs, with "and N more" where there are more.

### Changed

- The NAV page asks for the route's graph and draws it under the course diagram's route: a faint stub from each stop to each neighbour off the route (a few per stop, nearest first), "and N more" over a stop that has more, and a dotted arc, with its length, between two stops of one row that are neighbours but not a hop. Links between stops on different rows of a wrapped diagram are counted in a note.
- The panel and its legend say a link is a neighbour relation, not a travel lane. The route's own hops are always drawn over the neighbours, and each stub has a text alternative with its distance.

generation-output: unchanged
