**Headline:** Clicking a gas giant or habitable-world mark on the Galaxy Map opens its system instead of "no system with ID".

### Fixed
- **The Galaxy Map's marked systems (gas giant, habitable world and the other reasons) linked to the wrong ID.** The map nodes listed each marked system by its row id, but the system page and the nearest-systems list take the printed Object ID, so View system answered "no system with ID 345604". The nodes now carry the printed ID, and stored nodes are rebuilt once (node format 10).
