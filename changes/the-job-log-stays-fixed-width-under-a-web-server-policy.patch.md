**Headline:** The Generate log stays in a fixed-width font when the web server adds its own Content-Security-Policy.

### Fixed
- **The Generate log was no longer fixed width (UX.148).** Xterm.js sets the terminal rows' font in a `<style>` element it writes itself. A web server that sends its own Content-Security-Policy as well as the app's (Boss's Apache does) enforces both, and its policy, with no `'unsafe-inline'` for styles, blocked that element, so the log fell back to the page's proportional font and columns stopped lining up. `style.css` now gives `.job-terminal .xterm-rows` the monospace font, size and `white-space: pre` too. The deployment guide already says not to repeat the security headers in the web server; removing the server's own header is still the better setup. `test_the_log_stays_fixed_width_when_the_web_server_adds_a_stricter_policy` fails without the change.

TODO items resolved: UX.148. Added: none.
