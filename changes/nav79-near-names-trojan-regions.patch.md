**Headline:** The nearby-objects list no longer fails on a system that has a trojan region.

### Fixed
- Listing the bodies of a system for the nearby search raised a KeyError when the system had a giant's trojan region: the list read a belt number for every belt, and a trojan region takes none. It now names every belt the way the rest of the product does (`Asteroid belt N`, or `<host> L4 trojan region`).

generation-output: unchanged
generation-affects: none
