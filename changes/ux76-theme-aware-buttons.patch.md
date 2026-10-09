### Fixed
- A visited secondary link-button (outlined, such as Navigate or Show on Galaxy Map) kept the filled button's dark text on a clear background, so it was dark on dark in dark mode. Secondary links now keep their own colors visited or not.
- The outlined button's hover no longer drops its text below readable contrast in light mode; it thickens its outline instead of shading the background.
- A browser test checks every button look (filled, outlined, active, pressed, danger, map toggles) at rest, hovered and focused at 4.5:1 text contrast, and an icon from the sprite at 3:1, with the OS in dark or light and the site's own theme choice following or overriding it.
