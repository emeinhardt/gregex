# `gregex`

## Motivation / context

Glycans are a type of molecule that typically have a chain and/or tree-like structure. A number of researchers (e.g. Banin et al., 2002; Krambeck et al. 2009) have proposed compact, machine and (somewhat) human-readable notation ('linear code') for individual glycans, sets of formally similar glycans, and reactions. To represent formally similar sets of glycans, uncertainty operators (analogous to `.`, `+`, and `*` in regular expression tools) are employed. These operators do not currently have a mathematically precise or thorough explication, nor is there any open-source software that might fill a similar role.

This package contains functions (and a command-line interface to the main script) for clarifying and comparing the meaning of each of Krambeck et al. 2009's three uncertainty operators (`_`, `...`, `|`).

## Usage

 - The code in this repository can be imported as a package for programmatic use: `import gregex`.
 - The command-line interface can be accessed via the usual `python -m gregex ...` route. **`python -m gregex -h` will bring up the `argparse` help.**

### CLI example usage



## TODO

1. Migrate tests from the `dev` Jupyter notebook into `pytest` tests.
2. Add additional tests for code unique to `gregex.py` relative to the dev notebook.
3. Create a clean demo notebook from the existing development notebook.
4. Create a CLI for converting a linear code representation to an s-expression.

