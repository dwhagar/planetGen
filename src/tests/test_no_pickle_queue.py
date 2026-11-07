# tests/test_no_pickle_queue.py

"""
SEC.39: no queue carries pickle. RQ's default serializer is pickle, so a queue, a worker or a job fetched without a serializer would
read an object that runs code when a process that can write to Redis queues one. Every place that makes one names `JsonSerializer`
(`queue/redisqueue.py`). No database or Redis is needed.
"""

import ast
import os

from planetgen.queue import redisqueue

SOURCE = os.path.join(os.path.dirname(os.path.dirname(__file__)), "planetgen")


def _calls(tree, names):
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            function = node.func
            name = function.attr if isinstance(function, ast.Attribute) else getattr(function, "id", None)
            if name in names:
                yield node


def test_a_queue_made_by_the_site_serializes_with_json():
    for name in (*redisqueue.QUEUES, "planetgen-run-0123", "planetgen-web-0123"):
        assert redisqueue.queue(name, connection=object()).serializer is redisqueue.JsonSerializer


def test_nothing_that_makes_an_rq_queue_or_worker_leaves_the_serializer_to_rq():
    """A `Queue(...)`, `Worker(...)` or `Job.fetch(...)` without `serializer=` reads and writes pickle."""
    offenders = []
    for folder, _dirs, files in os.walk(SOURCE):
        for name in files:
            if not name.endswith(".py"):
                continue
            path = os.path.join(folder, name)
            tree = ast.parse(open(path, encoding="utf-8").read())
            for call in _calls(tree, {"TimedQueue", "WarmWorker", "RoundRobinWorker", "fetch"}):
                function = call.func
                if getattr(function, "attr", getattr(function, "id", None)) == "fetch" and not (
                        isinstance(function, ast.Attribute) and ast.unparse(function.value).endswith("Job")):
                    continue            # only RQ's `Job.fetch`, not another module's fetch
                if not any(keyword.arg == "serializer" for keyword in call.keywords):
                    offenders.append(f"{os.path.relpath(path, SOURCE)}:{call.lineno}")
    assert offenders == []


def test_the_queue_code_does_not_use_pickle():
    offenders = []
    for folder, _dirs, files in os.walk(os.path.join(SOURCE, "queue")):
        for name in files:
            if not name.endswith(".py"):
                continue
            for node in ast.walk(ast.parse(open(os.path.join(folder, name), encoding="utf-8").read())):
                imported = [alias.name for alias in node.names] if isinstance(node, ast.Import) else \
                    [node.module] if isinstance(node, ast.ImportFrom) else []
                if any(str(module).split(".")[0] in ("pickle", "cPickle", "dill", "cloudpickle", "marshal") for module in imported) \
                        or (isinstance(node, ast.Name) and node.id == "pickle"):
                    offenders.append(name)
    assert offenders == []
