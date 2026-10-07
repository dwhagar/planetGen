# stellarObjects/__init__.py

"""
Stellar Objects Package
=======================

The generation code still waiting to move into the `planetgen` package
(docs/design/library-migration.md, section 6). Its modules are imported
by name (`from planetgen.generation.star import Star`); the package itself
exports nothing, so importing one module never drags in the rest (which
would make import cycles with the moved `planetgen` modules).
"""
