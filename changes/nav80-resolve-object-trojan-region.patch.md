**Headline:** The page and the API answer for a giant's trojan region work, and a facility built on one names it.

### Fixed

- The object resolver (`/api/objects/belt:<id>`, the page of a belt and its parent chain) read a belt number for every belt, and a trojan region takes none, so it raised a KeyError for any region and its page failed. It names the region for its giant now (`<host> L4 trojan region`), as the lists of belts do. The same fix is made where a facility's host chain names a belt.

generation-output: unchanged
generation-affects: none
