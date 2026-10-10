### Changed
- The Galaxy Map and Sector pages load their scripts faster on a first visit (MAP.161): each page names the map's whole module tree in `modulepreload` hints, so the browser fetches the ~30 files in parallel instead of one after another, and the example Apache config explains how to turn on HTTP/2.
