**Headline:** A double click on the map is told by the time of its two clicks, so a first click the page is slow to answer no longer makes the second one a click again.

### Fixed
- Map controls: `mapcontrol.js` decided whether a click was the second of a pair by reading the clock when the page got to the event, with the first click's time taken before its handler ran. When the first click's handler (picking the star, drawing its panel) took longer than the 350 ms window, as it can on a busy machine, the second click came out as a single click again and the camera never flew; `test_a_double_click_on_a_star_flies_the_camera_to_it` failed that way once in a loaded gate (TEST.252: camera target unchanged after the double click). The pair is now told by the events' own time stamps. A new browser test stalls the first click for 500 ms and still sees the flight (it failed on the old code), and a unit test pins the rule.

TODO items resolved: TEST.252 / Added: none
