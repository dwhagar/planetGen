import os
import re

from setuptools import find_packages, setup

# NLTK's 'words' corpus (needed by planetgen.names.wordlists for name generation) is
# downloaded lazily, at import time, by planetgen/names/wordlists.py itself — it is
# NOT downloaded here. A previous version of this file tried to download it via
# a custom post-install command, but that imported nltk at the top of setup.py,
# which crashes: pip builds this package in an isolated build environment that
# only has setuptools available, not this package's own runtime dependencies
# (those are installed *after* the build, from install_requires below). That
# made `pip install .` fail unconditionally with ModuleNotFoundError: No module
# named 'nltk'.

here = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


def read_version():
    """
    Reads `__version__` out of `planetgen/_version.py` as plain text (a
    regex, not an import): the generator modules import `nltk` at module
    load time, and this script runs in an isolated build
    environment that doesn't have `nltk` (or any other runtime dependency)
    installed yet (see the note above on the corpus download).
    """
    version_path = os.path.join(here, 'src', 'planetgen', '_version.py')
    with open(version_path, encoding='utf-8') as f:
        contents = f.read()
    match = re.search(r"^__version__\s*=\s*['\"]([^'\"]+)['\"]", contents, re.M)
    if not match:
        raise RuntimeError("Unable to find __version__ in planetgen/_version.py")
    return match.group(1)


setup(
    name='planetGen',
    version=read_version(),
    # Everything lives under src/planetgen (src layout); src/tests and
    # src/html are not part of the package.
    packages=find_packages(where='src', include=['planetgen', 'planetgen.*']),
    package_dir={'': 'src'},
    package_data={
        # Data files the modules read at run time (setuptools leaves non-.py
        # files out of a package by default, and the module would crash on
        # first use): the common-password list (planetgen.admin.auth),
        # the offensive-word list (planetgen.names.wordlists) and the two
        # database schemas (planetgen.db.store).
        'planetgen.admin': ['common_passwords.txt.gz', 'common_passwords.LICENSE'],
        'planetgen.db': ['schema.sql', 'control_schema.sql', 'migrations/*.py', 'migrations/*.mako', 'migrations/*.md',
                         'migrations/versions/*.py'],
        'planetgen.names': ['offensive_words.txt'],
    },
    entry_points={
        'console_scripts': [
            'planetgen=planetgen.cli.generate:main',
        ],
    },
    install_requires=[
        # Lower-bounded (not pinned exact/upper-bounded) to floors with no
        # known security advisories as of this writing, so patch releases
        # keep flowing without needing this file edited for each one.
        'nltk>=3.9.1',
        # MySQL persistence (planetgen/db/store.py, TODO.md Phase 5) --
        # every entry point that touches the database (the `planetgen` command,
        # planetgen.db.query, planetgen.cli.migrate, the html/ CGI browser) needs these,
        # not just the Flask API, so they're core requirements rather
        # than the 'api' extra below.
        'pymysql>=1.1.1',
        # Admin password hashing (planetgen/admin/auth.py's lazy
        # `werkzeug.security` import) -- planetgen.cli.migrate calls
        # bootstrap_control_schema() unconditionally (install.sh/update.sh
        # both always run it), so this needs to be installed regardless of
        # whether the 'api' extra (the Flask API itself) is, same
        # reasoning as pymysql above. Floor matches the 'api'
        # extra's own flask>=3.0.3 (which already pulls in werkzeug>=3.0.0
        # transitively), so installing both extras together never
        # downgrades it.
        'werkzeug>=3.0.0',
        # Progress bar for `planetgen galaxy`'s sector-generation loops
        # (TimeElapsedColumn/TimeRemainingColumn; see _generation_progress
        # in planetgen/generation/run_common.py).
        'rich>=13.7.0',
        # The library migration (OPS.20, docs/design/library-migration.md):
        # pinned here first (OPS.21) so each swap PR only changes code.
        # The work queue and web jobs on Redis (PERF.24), and Redis as
        # Flask-Limiter's storage (SEC.30).
        'redis>=5.0.0',
        'rq>=1.16.0',
        # Two-step sign-in and its QR code (SEC.29).
        'pyotp>=2.9.0',
        'segno>=1.6.0',
        # Markdown rendering (UX.39).
        'markdown>=3.6',
        # The page and tile caches (PERF.25). diskcache is left out: its
        # last release has an unfixed advisory (PYSEC-2026-2447, unpickling
        # whatever is in the cache folder), which the dependency audit
        # rejects; PERF.25 decides what replaces tilecache.py.
        'cachetools>=5.3.0',
        # The database layer and migrations (DB.11).
        'sqlalchemy>=2.0.30',
        'alembic>=1.13.0',
        # Input validation (ADM.21).
        'pydantic>=2.7.0',
        # Physics (GEN.66) and nebula meshes (marching cubes).
        'numpy>=1.26.0',
        'scipy>=1.13.0',
        'astropy>=6.0.0',
        'scikit-image>=0.22.0',
    ],
    extras_require={
        # matplotlib/numpy are for tests/galaxy_shape_visualizer_cli.py
        # alone (renders the density model as an actual image) -- no
        # other test or application code touches either.
        'test': ['pytest>=7.4.0', 'matplotlib>=3.8.0', 'numpy>=1.26.0',
                 # tests/test_fuzz_*.py (property-based brute-force tests;
                 # see tests/fuzz_support.py).
                 'hypothesis>=6.100.0',
                 # `pytest -n auto` (docs/testing.md); the suite still
                 # runs serially without it.
                 'pytest-xdist>=3.5.0',
                 'pytest-timeout>=2.3.1'],
        'api': ['flask>=3.0.3', 'flask-limiter>=3.7.0'],
        # tests/test_web_a11y.py alone (headless-browser layout and
        # accessibility checks of the Flask pages); it skips without it.
        # Chromium itself comes from `python -m playwright install chromium`.
        'browser': ['playwright>=1.49.0'],
        # The standalone app server where Apache's mod_wsgi isn't used:
        # gunicorn under launchd on macOS (docs/deployment/macos.md).
        # install.sh on macOS installs it from requirements-server.lock.
        'server': ['gunicorn>=23.0.0'],
    },
    author='David Hagar',
    author_email='david.hagar@gmail.com',
    description='A procedural planet and star system generator.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    license='CC0-1.0',
    license_files=('LICENSE.md',),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: CC0 1.0 Universal (CC0 1.0) Public Domain Dedication',
        'Operating System :: OS Independent',
    ],
    # Driven by the 'api' extra's floors above (Flask 3.x/PyMySQL 1.1.x
    # both need 3.9+) rather than anything in this package's own code.
    python_requires='>=3.9',
)