# tests/test_no_pickle_queue.py

"""
SEC.39: no queue carries pickle: a process that can write to Redis must not be able to queue an object that runs code when it is read. OPS.235: RQ (whose default serializer was pickle) is gone, so
what is left to check is that nothing under `queue/` imports a pickling module. No database or Redis is needed.
"""

import ast
import os

SOURCE = os.path.join(os.path.dirname(os.path.dirname(__file__)), "planetgen")


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
