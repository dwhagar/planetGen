**Headline:** Bringing a galaxy to a new epoch no longer generates the empty cells rogue planets landed in.

### Changed
- **bring-to-epoch leaves uncharted cells alone (GEN.1354).** A drifting rogue planet that lands in a cell nothing was generated in leaves a placeholder sector row there. The run took those rows for sectors made under an older epoch and generated each of them: in a small test galaxy 74 of the 75 sectors it would have redone were placeholders, and every one generated charted a whole new sector whose own rogues landed in more empty cells. The run now redoes only the sectors that were charted, counts only those as the fill stage when it reports what is behind, and says how many uncharted cells it left alone. A placeholder keeps its rogues and is generated when someone charts the cell.

generation-output: unchanged

TODO items resolved: GEN.1354.
