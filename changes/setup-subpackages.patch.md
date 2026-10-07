### Fixed
- **`pip install` failed after the first package move.** `setup.py` mapped `planetgen.util` to a folder named `src/planetgen.util`; subpackages now map to their real folders.
