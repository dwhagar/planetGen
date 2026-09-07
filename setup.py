import os
import re

from setuptools import find_packages, setup

# NLTK's 'words' corpus (needed by stellarObjects.names for name generation) is
# downloaded lazily, at import time, by stellarObjects/names.py itself — it is
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
    Reads `__version__` out of `stellarObjects/_version.py` as plain text
    (regex, not import) rather than `import stellarObjects` -- that package's
    `__init__.py` imports `Planet`, which imports `names.py`, which imports
    `nltk` at module load time, and this script runs in an isolated build
    environment that doesn't have `nltk` (or any other runtime dependency)
    installed yet (see the note above on the corpus download).
    """
    version_path = os.path.join(here, 'stellarObjects', '_version.py')
    with open(version_path, encoding='utf-8') as f:
        contents = f.read()
    match = re.search(r"^__version__\s*=\s*['\"]([^'\"]+)['\"]", contents, re.M)
    if not match:
        raise RuntimeError("Unable to find __version__ in stellarObjects/_version.py")
    return match.group(1)


setup(
    name='planetGen',
    version=read_version(),
    packages=find_packages(),
    package_data={
        # stellarObjects.names reads this at import time; setuptools does not
        # include non-.py files in a package by default, so without this the
        # installed package is missing the file and crashes on first import.
        'stellarObjects': ['offensive_words.txt'],
    },
    py_modules=['systemGen', 'sectorGen'],
    entry_points={
        'console_scripts': [
            'systemgen=systemGen:main',
            'sectorgen=sectorGen:main',
        ],
    },
    install_requires=[
        'nltk',
    ],
    extras_require={
        'test': ['pytest'],
        'api': ['flask'],
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
    python_requires='>=3.6',
)