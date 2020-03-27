# `gregex`
A tool for investigating and working with regular-expression-like operators that describe glycans.

## Motivation / context

Glycans are a type of molecule that typically have a chain and/or tree-like structure. A number of researchers (e.g. Banin et al., 2002; Krambeck et al. 2009) have proposed compact, machine and (somewhat) human-readable notation ('linear code') for individual glycans, sets of formally similar glycans, and reactions. To represent formally similar sets of glycans, uncertainty operators (analogous to `.`, `+`, and `*` in regular expression tools) are employed. These operators do not currently have a mathematically precise or thorough explication, nor is there any open-source software that might fill a similar role.

This package contains functions (and a command-line interface to the main script) for clarifying and comparing the meaning of each of Krambeck et al. 2009's three uncertainty operators (`_`, `...`, `|`) as documented there and in Glymmer manual.

## Usage

 - The code in this repository can be imported as a package for programmatic use: `import gregex`.
 - The command-line interface can be accessed via the usual `python -m gregex ...` route. **`python -m gregex -h` will bring up the `argparse` help.**
    - Only a fraction of the package's functionality is currently exposed through the command-line interface.

### CLI example usage

Provided the `gregex` module is on your path (via e.g. step 2 of the installation process below), some of the functionality of the `gregex` Python module is available via `gregex -m gregex <ARGS>`. 

All CLI functionality performs some operation on a single glycan linear code expression (the first and main argument to the script). Exactly which operation is dictated by other flags and arguments.

For complete details and a description of *all* functionality and flags, use the command-line help flag: `python -m gregex -h`.

#### Converting a linear code representation of a glycan to an s-expression

The [native representation](https://en.wikipedia.org/wiki/S-expression) of code and data in Lisp dialects makes tree and list structure readily apparent, even for longer glycans, particularly when indented according to common conventions.

`python -m gregex 'NNa6Ab4GNb4(NNa3(ANb4)Ab4GNb2)Ma3(NNa3(ANb4)Ab4GNb3Ab4GNb2(NNa3(ANb4)Ab4GNb6)Ma6)Ma4GNb4(Fa6)GN' -e`

(currently) yields

`(GN Fa6 (GNb4 (Ma4 (Ma6 (GNb6 (Ab4 ANb4 NNa3)) (GNb2 (Ab4 (GNb3 (Ab4 ANb4 NNa3))))) (Ma3 (GNb2 (Ab4 ANb4 NNa3)) (GNb4 (Ab4 NNa6))))))`

(`gregex` currently doesn't do pretty-printing of s-expressions, but for the time being, any widely-used text editor will support packages that automatically indent s-expressions according to common conventions.)

#### Checking whether a string matches an uncertainty operator in some glycan

`python -m gregex 'Ma6_M' -s '(Ma4)'` checks whether `(Ma4)` can be substituted for the operator `_` in `Ma6_M`. It can, so this returns `True` to stdout.

#### Finding all nonempty matches for your choice of uncertainty operator in some glycan

`python -m gregex 'Ma6(Ma4)M' -o '_'` writes a set of lines to stdout indicating all the nonempty subsequences of `Ma6(Ma4)M` that could be replaced with `_` and yield a syntactically well-formed linear code expression.

`python -m gregex 'Ma6(Ma4)M' -o '_' -c` is the same, but each line now contains
	left_context	match	right_context
for some match. 

`python -m gregex 'Ma6(Ma4)M' -o '_' -s '(Ma2)' -c` is similar to the previous command, but checks for each `(left context, match, right_context)` triple whether `(Ma2)` can successfully match the location of `_` in each possible left-match-right split of the original linear code expression.

## Requirements / installation

All code has been developed and tested on Ubuntu 18.04.3 and MacOS 10.13.5.

The three most salient dependencies are
 - [`funcy`](https://funcy.readthedocs.io/en/stable/)
 - [`glypy`](https://pythonhosted.org/glypy/) (so far only necessary for development, not for CLI functionality or most other functions)
 - `Python 2.7`
    - `glypy` does not currently support Python 3.

To set up a new conda environment that contains this repository's dependencies,
1. `git clone` this repository to a filepath of your choice.
2. `cd path_to_repo`
3. Create the conda environment automatically via the `.yml` file in the repository (`conda env create -f gregex_env.yml`, followed by `conda activate gregex`) *or* enter the commands in `conda_manual_environment_creation.txt` at your command prompt, one at a time.

## TODO

1. Migrate tests from the `dev` Jupyter notebook into `pytest` tests.
2. Add additional tests for code unique to `gregex.py` relative to the dev notebook.
3. Create a clean demo notebook from the existing development notebook.
4. Add pretty-printing support to s-expression conversion and make argument labels (=bond information) more explicit. For example: 
   
`(GN Fa6
    (GNb4 (Ma4 (Ma6 (GNb6 (Ab4 ANb4
                               NNa3))
                    (GNb2 (Ab4 (GNb3 (Ab4 ANb4
                                          NNa3)))))
               (Ma3 (GNb2 (Ab4 ANb4
                               NNa3))
                    (GNb4 (Ab4 NNa6))))))`

`(GN :a6 F
    :b4 (GN :a4 (M :a6 (M :b6 (GN :b4 (A :b4 AN
                                         :a3 NN))
                          :b4 (GN :b4 (A :b3 (GN :b4 (A :b4 AN
                                                        :a3 NN)))))
                   :a3 (M :b2 (GN :b4 (A :b4 AN
                                         :a3 NN))
                          :b4 (GN :b4 (A :a6 NN))))))`

[`hy`](https://docs.hylang.org/en/stable/) plausibly has pretty-printing facilities that support this out-of-the-box.

