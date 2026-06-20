"""
Shared SQLite sink for all Final_plots figure scripts.

Each figure writes the exact data it plots into a table in plot_data.db so the
underlying numbers are archived and queryable alongside the PNGs.
"""
from pathlib import Path
import sqlite3
import pandas as pd

DB_PATH = Path(__file__).parent / 'plot_data.db'


def write_table(df: pd.DataFrame, table: str, db_path: Path = DB_PATH):
    """Replace `table` in the figures DB with the contents of `df`."""
    df = df.copy()
    with sqlite3.connect(db_path) as conn:
        df.to_sql(table, conn, if_exists='replace', index=False)
    print(f'  [db] wrote {len(df)} rows -> {table} ({db_path.name})')
