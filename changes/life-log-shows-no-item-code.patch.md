**Headline:** The sealed-ocean oxygen line of the generation log no longer names a TODO item.

### Fixed
- `generation/life.py`: the log choice "Sealed ocean oxygen" ended its reason with "(GEN.1824)", a TODO code in a string the code shows, so `test_no_todo_codes_shown` was red on main after Train 54. The code is dropped; the reason reads the same otherwise.

generation-output: unchanged
