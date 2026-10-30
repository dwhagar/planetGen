**Headline:** Pressing Back twice quickly across a sector change now lands on the right system.

### Fixed
- **Two quick presses of Back across a sector change (MAP.352).** Pressing Back twice in quick succession, the second time to a system in a different sector than the one the first press was moving to, lost the system: the second entry's open request was handed to the display the first move was arriving at and dropped. A move still under way is now finished before the next history entry is read, and a sector or system still loading for the view being left is cancelled when the next view is taken up, so the second press lands on its system.
