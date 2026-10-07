# planetgen/cli/worker.py

"""
One RQ worker for planetGen's queues (PERF.24), at the lowered priority
generation workers always run at. `planetgen.queue.redisqueue.
ensure_workers` starts these in burst mode, so each exits once its
queues are empty.

Usage:
    python -m planetgen.cli.worker [--burst] [--name NAME] [--url REDIS_URL] [QUEUE ...]
"""

import argparse

from planetgen.queue import redisqueue, work


def main(argv=None):
    parser = argparse.ArgumentParser(prog="planetgen.cli.worker", description=__doc__.strip().splitlines()[0])
    parser.add_argument("queues", nargs="*", default=list(redisqueue.QUEUES),
                        help="the queues to serve, first one first (default: all of planetGen's)")
    parser.add_argument("--burst", action="store_true", help="exit once the queues are empty")
    parser.add_argument("--name", help="the worker's name (default: RQ picks one)")
    parser.add_argument("--url", help="the Redis server (default: PLANETGEN_REDIS_URL, else config.json's redis.url)")
    args = parser.parse_args(argv)
    work.lower_priority()
    connection = redisqueue.connect(args.url)
    queues = [redisqueue.queue(name, connection) for name in args.queues]
    worker = redisqueue.worker_class()(queues, connection=connection, name=args.name)
    worker.work(burst=args.burst, with_scheduler=False)


if __name__ == "__main__":
    main()
