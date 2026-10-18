**Headline:** The Stop waiting button on a page that is being built now works from the moment it appears.

### Fixed
- `static/slowload.js` awaited two dynamic imports (the number and duration formatters, loaded by the script's own version) at the top of the module, so its click handlers were attached only after both had arrived, while the page already showed the button. A click in that gap did nothing: the status line kept its estimate and the page went on polling. The buttons are now wired at once, the formatters load alongside, and the polling starts when they have arrived (Stop waiting, Check again and Try again act on the first click). Found as `test_stop_waiting_stops_the_polling_and_check_again_picks_it_up[dark-phone]` and `[light-phone]` failing under load while passing alone; the click event reached the button in every failing run (recorded), and the handler had not been attached. With the fix the test passed 80 runs in 80 under 12 busy loops (about one invocation in five failed before).

TODO items resolved: TEST.248 / Added: none
