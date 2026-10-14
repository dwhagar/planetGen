### Fixed
- Linking a big run's sectors to their neighbours no longer looks frozen: it works in batches of 20 neighbouring sectors instead of 200 scattered ones (a batch was about ten minutes on a dense run, with the progress bar standing still and the neighbour lock held throughout), so the bar moves every minute or so and the lock is short. The links stored are the same.
