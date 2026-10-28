**Headline:** The admin statistics and queue tiles now say what they count, and every numeric setting names its unit.

### Changed
- **Admin tiles and settings units (UX.212).** The Stats page tiles read "12 sectors", "3,456 systems", "2 names" instead of a bare number, and the Queue page's "Workers active" and "Runs live" tiles read "4 workers" and "1 run". The nine numeric settings that had no unit (database port, the three proxy counts, the two password hashing figures, cached answers, finished jobs kept, settings file version) now carry one. Tests pin both: no tile is a bare number, and every numeric setting has a unit.
