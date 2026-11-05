**Headline:** The settings review-page width test waits for the page's components before it measures the page.

### Fixed
- **A phone-width overflow reported once (606 px at the gear dropdown) came from measuring before the components had loaded.** The settings browser tests now wait until the dropdown components are defined (at most five seconds) before reading the page width, and a test holds the components script back two seconds to pin it.
