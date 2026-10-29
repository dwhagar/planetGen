**Headline:** A sector's and a system's scene now come up because the camera is near them, in any view, not only in a window of single sectors.

### Changed
- **Scenes follow the camera (MAP.335).** The camera flown close to a star from the whole galaxy, or from a window of bigger blocks, now brings up the sector's stars and then the star's system out of where the eye is, with no window of sectors in front of them. The map asks the server for the generated sectors near the eye (a new read, `/galaxy/near`, cut from the sector files the map already caches) and the wheel opens the sector and the system where the camera is.
