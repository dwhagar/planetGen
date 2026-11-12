**Headline:** The one-off system page can ask for a planet class or a life class and shows what came of it.

### Changed
- **Ask for a class and Ask for a life class on the one-off system page (UX.302).** Two boxes under Options send `--force-class` and `--force-life-class` (GEN.1620). The result shows a line for each request above the system: which body holds the class and how well its conditions fit, or that none came out. A bad letter or life class, or a class letter together with a system file's slots, is refused on the form, and a rogue planet or comet result refuses them as it does the other star-only options. The check of a class request moved from the command line to `generation/run_system.py` so the page and the command share it.

generation-output: unchanged
