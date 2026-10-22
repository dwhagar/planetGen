**Headline:** A zoom on the Galaxy Map no longer stalls the page walking every star for every tile that has not arrived.

### Changed
- **Carrying stars over to tiles still on the way (MAP.249).** While the tiles of a new zoom level arrive, the map keeps drawing the stars already on screen where the new tiles will be. For that it walked all the drawn stars once for each missing tile: about 250 tiles by 45,000 stars is 11 million box tests, 665 ms of main-thread time on the test machine, repeated on each fetch step of a zoom. It now places each star in its tile with arithmetic in one pass (6.5 ms for the same case) and carries the same stars.
- **A view known from its tiles (MAP.249).** Each time tiles arrived or the camera settled with every tile in hand, the map built a list of every star in the view (up to 45,000) and sorted and joined it, about 160 ms of main-thread time on each call, only to find the stars were the ones already drawn. It now compares the tiles themselves (`static/tileview.js`): a view whose tiles are the ones already drawn costs a pass over about 250 tile keys. Nothing a visitor sees changes; the page works the same with fewer pauses on a zoom.
- **Guards.** A test fails if a view asks for more than 220 tiles at any of eight camera distances from 30 to 100,000 parsecs, and a script test fails if signing a view of 250 tiles reads a star or takes over 5 ms.

TODO items resolved / Added: none / Added: none
