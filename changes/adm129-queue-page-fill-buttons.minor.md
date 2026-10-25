**Headline:** The Queue page can start the habitability re-score and the map-luminosity fill with a button.

### Added
- **Start a fill panel on the admin Queue page (ADM.129).** "Regenerate habitability" re-scores every planet and moon from
  its stored conditions (it re-draws nothing), and "Fill map luminosity" fills the map luminosity of scattered
  phenomena that lack it. Each button asks for confirmation, then starts `planetgen rescore-habitability` or
  `planetgen fill-map-luminosity` as a job with the rules of any generation run: refused while another run is live or
  everything is stopped, the one-run lock, an activity-log entry, progress on the Queue page, Stop everything and Retry.
