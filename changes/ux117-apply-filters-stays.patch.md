**Headline:** The Apply filters button no longer disappears from filtered tables, and it works whether or not the table has loaded.

### Changed
- **Apply filters stays on every filtered table.** Once the table had fetched its first rows the button was hidden, so a table whose load was slow or had failed gave no way to apply the filters you had ticked. The button now stays visible on the phenomena, sector contents and API keys tables.
- **It submits to the server until the table is live.** While the first rows are loading, or if that load fails, the button sends the filters to the server and shows the filtered page as it always did. Once the table is live it applies the filters in place, and if the last fetch failed it asks for the rows again.

TODO items resolved: UX.117 / Added: none
