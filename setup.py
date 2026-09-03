import os

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

setup(
    name='planetGen',
    version='0.1.0',
    packages=find_packages(),
    package_data={
        # stellarObjects.names reads this at import time; setuptools does not
        # include non-.py files in a package by default, so without this the
        # installed package is missing the file and crashes on first import.
        'stellarObjects': ['offensive_words.txt'],
    },
    py_modules=['planetGen'],
    entry_points={
        'console_scripts': [
            'planetgen=planetGen:main',
        ],
    },
    install_requires=[
        'nltk',
    ],
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