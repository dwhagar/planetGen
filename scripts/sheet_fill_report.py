#!/usr/bin/env python3
# scripts/sheet_fill_report.py

"""
The fill audit of the property sheet (UX.208): generates sample systems into a scratch database, counts for every kind of
body, class and sheet column how often it is filled and how often empty, and says for each empty one whether
`describe.sheet.APPLIES` says the column applies there. Prints the table and writes it to `reports/`.

    python scripts/sheet_fill_report.py [--systems 150] [--seed 1000] [--out /mnt/project-files/reports/sheet-fill-report.md]
"""

import argparse
import pathlib
import sys
import uuid

import pymysql

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]

from planetgen.db import store  # noqa: E402
from planetgen.db.store import MySQLConfig  # noqa: E402
from tests import sheet_fill_support as support  # noqa: E402


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--systems", type=int, default=150, help="how many systems to generate")
    parser.add_argument("--seed", type=int, default=1000, help="the run seed of the first system")
    parser.add_argument("--out", default=str(ROOT / "reports" / "sheet-fill-report.md"), help="where the table is written")
    args = parser.parse_args(argv)

    server = support.server_kwargs()
    name = f"planetgen_fill_{uuid.uuid4().hex[:12]}"
    admin = pymysql.connect(**server)
    with admin.cursor() as cur:
        cur.execute(f"CREATE DATABASE `{name}`")
    admin.commit()
    config = MySQLConfig(database=name, **server)
    try:
        store.get_connection(config).close()    # applies the schema
        support.generate(config, args.systems, args.seed)
        text = "\n".join(support.lines(support.tally(config)))
    finally:
        store.close_pool(config)
        with admin.cursor() as cur:
            cur.execute(f"DROP DATABASE `{name}`")
        admin.commit()
        admin.close()
    print(text)
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(f"# Property sheet fill report ({args.systems} systems from seed {args.seed})\n\n```\n{text}\n```\n",
                   encoding="utf-8")


if __name__ == "__main__":
    main()
