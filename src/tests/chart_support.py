# tests/chart_support.py

"""Helper for fixtures that insert `sectors` rows by hand: a saved sector is charted (`store.save_sector` records its cell in the chart ledger), and a placed row with no ledger row is an uncharted placeholder (MAP.399, MAP.402) that no block, region or count includes."""


def chart_all(conn):
    """Records the chart ledger row of every placed sector row that has none, as a saved sector would have."""
    conn.execute("INSERT IGNORE INTO chart_ledger (shell_index, band_index, wedge_index, generator_epoch)"
                 " SELECT shell_index, band_index, wedge_index, 1 FROM sectors WHERE shell_index IS NOT NULL")
