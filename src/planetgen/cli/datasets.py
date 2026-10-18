#!/usr/bin/env python3
# planetgen.cli.datasets

"""
Downloads the static data the generator reads (GEN.548), or says where each dataset stands.

    python3 -m planetgen.cli.datasets fetch [KEY ...]
    python3 -m planetgen.cli.datasets status [KEY ...]

`fetch` is the setup step `update.sh` and `install.sh` run as root: each file is downloaded once with the standard
User-Agent, checked against the catalogue's sha256 and kept in the cache folder (`PLANETGEN_DATA_DIR`, else
`/var/cache/planetgen/datasets`). It prints one line per dataset and exits 1 when any is still missing, which
the scripts report as a warning, never as a failure: a host that is down cannot stop an update. `status` fetches
nothing and exits 1 the same way: it is the data check in front of every Generate job (GEN.590). Nothing in the generator or the site downloads at run time.
"""

import argparse
import sys

from planetgen.datasets import catalogue, fetch


_CONSEQUENCE = " -- generation will not start until it is fetched (it replaces: {})"


def _chosen(keys):
    by_key = {dataset.key: dataset for dataset in catalogue.ALL}
    unknown = [key for key in keys if key not in by_key]
    if unknown:
        raise SystemExit(f"error: unknown dataset {', '.join(unknown)} (known: {', '.join(by_key) or 'none'})")
    return [by_key[key] for key in keys] or list(catalogue.ALL)


def _line(dataset, missing, text):
    return f"{dataset.key}: {text}" if not missing else f"{dataset.key}: missing ({dataset.fill}), {text}"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("command", choices=("fetch", "status"))
    parser.add_argument("keys", nargs="*", metavar="KEY", help="Datasets to handle (default: all)")
    args = parser.parse_args(argv)
    datasets = _chosen(args.keys)
    if not datasets:
        print("No datasets are registered yet.")
        return 0
    absent = 0
    for dataset in datasets:
        if args.command == "fetch":
            result = fetch.fetch(dataset)
            if result.status == "missing":
                absent += 1
                print(_line(dataset, True, "; ".join(result.reasons)) + _CONSEQUENCE.format(dataset.replaces))
            else:
                print(_line(dataset, False, f"{result.status} ({len(dataset.pieces)} files, {fetch.folder(dataset)})"))
        else:
            located = fetch.locate(dataset)
            if located.missing:
                absent += 1
                print(_line(dataset, True, located.reason) + _CONSEQUENCE.format(dataset.replaces))
            else:
                print(_line(dataset, False, f"present in {located.folder}"))
    return 1 if absent else 0


if __name__ == "__main__":
    sys.exit(main())
