#!/usr/bin/env python3
# planetgen.cli.datasets

"""
Downloads the static data the generator reads (GEN.548), or says where each dataset stands.

    python3 -m planetgen.cli.datasets fetch [KEY ...]
    python3 -m planetgen.cli.datasets status [KEY ...]

`fetch` is the setup step `update.sh` and `install.sh` run as root: each file is downloaded once with the standard
User-Agent, checked against the catalogue's sha256 and kept in the cache folder (`PLANETGEN_DATA_DIR`, else
`/var/cache/planetgen/datasets`). It prints one line per dataset, and under it one line per file with its https URL
and what happened to it (`fetched`, `kept` when it was already there and not requested again, or `failed` with the
reason), and exits 1 when any dataset is still missing, which the scripts report as a warning, never as a failure: a host
that is down cannot stop an update. `status` fetches nothing and exits 1 the same way: it is the data check in front of
every Generate job (GEN.590), and it prints the URL of every file that is missing. Nothing in the generator or the site
downloads at run time.
"""

import argparse
import sys

from rich.console import Console
from rich.progress import BarColumn, DownloadColumn, Progress, TextColumn, TimeRemainingColumn, TransferSpeedColumn

from planetgen.datasets import catalogue, fetch
from planetgen.util import log


_CONSEQUENCE = " -- generation will not start until it is fetched (it replaces: {})"
_DISREGARDED = " -- the run takes a labelled fallback: {}"


def _chosen(keys):
    by_key = {dataset.key: dataset for dataset in catalogue.ALL}
    unknown = [key for key in keys if key not in by_key]
    if unknown:
        raise SystemExit(f"error: unknown dataset {', '.join(unknown)} (known: {', '.join(by_key) or 'none'})")
    return [by_key[key] for key in keys] or list(catalogue.ALL)


def _line(dataset, missing, text):
    return f"{dataset.key}: {text}" if not missing else f"{dataset.key}: missing ({dataset.fill}), {text}"


def _piece_line(state, piece, reason=""):
    reason = reason.removeprefix(piece.url + ": ")      # the download errors already start with the URL
    return f"    {state:<8} {piece.url}" + (f" ({reason})" if reason else "")


class _Bars:
    """
    One bar per file while it downloads (the name, the bytes done of the total, the speed and the time left, in the style of the generation
    bar of `generation.run_common`), taken away when the file is done; the `fetch` line below it says how it went. Lines printed while a
    bar is live go above it (`note`).
    """

    def __init__(self, console):
        self._progress = Progress(
            TextColumn("[progress.description]{task.description}"), BarColumn(), DownloadColumn(), TransferSpeedColumn(),
            TextColumn("[dim]remaining"), TimeRemainingColumn(), console=console, transient=True)
        self._tasks = {}

    def begin(self, piece):
        self._progress.start()
        self._tasks[piece.name] = self._progress.add_task(piece.name, total=piece.size)

    def update(self, piece, count):
        self._progress.update(self._tasks[piece.name], completed=count)

    def note(self, text):
        self._progress.console.print(f"    {text}", markup=False, highlight=False)

    def end(self, piece):
        task = self._tasks.pop(piece.name, None)
        if task is not None:
            self._progress.remove_task(task)
        if not self._tasks:
            self._progress.stop()


def _download_progress():
    """The bars on a terminal, plain periodic lines anywhere else (a log file, a pipe)."""
    console = Console()
    if console.is_terminal:
        return _Bars(console)
    return fetch.PlainProgress(lambda text: print(text, flush=True))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("command", choices=("fetch", "status"))
    parser.add_argument("keys", nargs="*", metavar="KEY", help="Datasets to handle (default: all)")
    parser.add_argument("--disregard-missing-data", action="store_true",
                        help="status: say that a missing dataset with a fallback is taken as one and exit 0 for it (the "
                             "Generate page's data check when the run was told to disregard missing data)")
    log.add_debug_arguments(parser)
    args = parser.parse_args(argv)
    log.apply_debug_arguments(args, parser)
    datasets = _chosen(args.keys)
    if not datasets:
        print("No datasets are registered yet.")
        return 0
    absent = 0
    progress = _download_progress() if args.command == "fetch" else None
    for dataset in datasets:
        if args.command == "fetch":
            result = fetch.fetch(dataset, progress=progress)
            if result.status == "missing":
                absent += 1
                print(_line(dataset, True, "; ".join(result.reasons)) + _CONSEQUENCE.format(dataset.replaces))
            else:
                print(_line(dataset, False, f"{result.status} ({len(dataset.pieces)} files, {fetch.folder(dataset)})"))
            for handled in result.pieces:
                print(_piece_line(handled.state, handled.piece, handled.reason))
        else:
            located = fetch.locate(dataset)
            if located.missing:
                taken = args.disregard_missing_data and bool(dataset.fallback)
                absent += not taken
                print(_line(dataset, True, located.reason)
                      + (_DISREGARDED.format(dataset.fallback) if taken else _CONSEQUENCE.format(dataset.replaces)))
                for piece, reason in located.lacking:
                    print(_piece_line("missing", piece, reason))
            else:
                print(_line(dataset, False, f"present in {located.folder}"))
    return 1 if absent else 0


if __name__ == "__main__":
    sys.exit(main())
