# `gregex`
A tool for investigating and working with regular-expression-like operators that describe glycans in linear code. Currently an alpha release.

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

#### Checking syntactic well-formedness of linear code representations

`python -m gregex 'Ma6(Ma4)M'` returns a boolean indicating the linear code expression is well-formed or not according to the following [context-free grammar](https://en.wikipedia.org/wiki/Context-free_grammar): 

```
exp ⟶ subexp non_main_branch+ stem | stem | λ
stem ⟶ SU_with_bond_info* SU_bare
non_main_branch ⟶ '(' subexp ')'
subexp ⟶ substem | subexp non_main_branch+ substem
substem ⟶ SU_with_bond_info+
SU_with_bond_info ⟶ SU_bare bond_type bond_location
bond_type ⟶ 'a' | 'b' | '?'
bond_location ⟶ '1' | '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9' | '?'
SU_bare ⟶ 'A' | 'AN' | 'B' | 'E' | 'F' | 'G' | 'GN' | 'G[Q]' | 'H' | 'H[2Q, 4Q]' | 'I' | 'K' | 'L' | 'M' | 'NG' | 'NJ' | 'NN' | 'NN[9N]' | 'N[5Q]' | 'O' | 'P' | 'PH' | 'R' | 'S' | 'U' | 'W' | 'X'
```


where
 - `⟶`, `|`, `λ` , `*`, and `+` are all reserved and/or metalinguistic symbols with their usual formal-language theoretic meaning (see any textbook or introductory material for reference).
 - all terminal symbols are quoted string literals, except for the empty string.
 - the enumeration of saccharide units is taken from a relatively arbitrary mix of what `glypy` supports and what `glymmer` supports.

Note that this is a declarative specification of what linear code expressions are that represent a single glycan or a set of glycans (via the uncertainty operators about bond type and position). See e.g. [BNF](https://en.wikipedia.org/wiki/Backus%E2%80%93Naur_form) for more on why specifications like this are common.

**NOTE 1:** The parser does *exactly* what it says on the tin: it checks syntactic well-formedness. It is somewhere between *not* the job of a parser (or not something you really want a parser to do) to check or enforce things like:
 - syntactic *conventions* about the linear ordering of children
 - whether a linear code representation describes something physically possible (= the *denonational semantics* of linear code).

Some other part of `gregex` might support these features on top of parsing eventually, but for now they are absent. (`glypy` might support some aspects of the second feature.) 

**NOTE 2:** As you may have noticed, with the exception of bond type/location uncertainty operators, linear code expressions with uncertainty operators are *not* part of this grammar. This is a consequence of their current ad-hoc definition in terms of string-matching. Incorporating them into the parser is possible through ad-hoc hacks and further research clarifying their meaning. 

**NOTE 3:** For longer linear code expressions (e.g. the large glycan example elsewhere on this page), current code may need a few tens of GB and a few minutes to calculate well-formedness. While this is not a problem for servers commonly used in scientific computing, it may not be practical for use on a researcher's personal laptop. Since NLTK is largely a research and pedagogically-oriented library, a more performant parser could easily improve on this.

#### Converting a linear code representation of a glycan to an s-expression

While linear code is more compact than more general tree notations when chaining ('unary branching') is more typical than (multi-child) branching, the 'bushier' a glycan is and the more monosaccharides are in the glycan, the harder it will be for a human to see hierarchical structure at a glance and the more likely they are to make mistakes while reading or editing. 
`gregex` has (currently somewhat limited) support for exporting a glycan represented in linear code to a notation that makes the tree structure more apparent: 's-expressions'.

[This representation](https://en.wikipedia.org/wiki/S-expression) of code and data native to Lisp dialects makes tree and list structure readily apparent, even for longer glycans, particularly when indented according to common conventions.
S-expressions ('s-exps') also have a long history of use in natural language parsing for creating human- and machine-readable representations of syntactic trees.

```
python -m gregex 'NNa3(ANb4)Ab4GNb2(NNa6Ab4GNb4)Ma3(NNa3(ANb4)Ab4GNb3Ab4GNb2(NNa3(ANb4)Ab4GNb6)Ma6)Ma4GNb4(Fa6)GN' -e
```

(currently) yields

```
(GN Fa6 (GNb4 (Ma4 (Ma6 (GNb6 (Ab4 ANb4 NNa3)) (GNb2 (Ab4 (GNb3 (Ab4 ANb4 NNa3))))) (Ma3 (GNb4 (Ab4 NNa6)) (GNb2 (Ab4 ANb4 NNa3))))))
```

(`gregex` currently doesn't do pretty-printing of s-expressions, but for the time being, any widely-used text editor will support packages that automatically indent s-expressions according to common conventions. See the `TODO` item below for how this pretty-printed output would likely appear.)

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

The four most salient dependencies are
 - [`funcy`](https://funcy.readthedocs.io/en/stable/), supporting functional programming.
 - [`nltk`](https://www.nltk.org/), for linear code expression parsing outside of `glypy`.
 - [`glypy`](https://pythonhosted.org/glypy/) (so far only necessary for development, not for CLI functionality or most other functions)
 - `Python 2.7`
    - `glypy` does not currently support Python 3.
        - `gregex` should otherwise be Python 3 compatible.
    - Note that nearly every direction for further development of this package depends on third-party packages that have at best limited support for Python 2.

To set up a new conda environment that contains this repository's dependencies,
1. `git clone` this repository to a filepath of your choice.
2. `cd path_to_repo`
3. Create the conda environment automatically via the `.yml` file in the repository (`conda env create -f gregex_env.yml`, followed by `conda activate gregex`) *or* enter the commands in `conda_manual_environment_creation.txt` at your command prompt, one at a time.

### Optional/complementary packages

[`csvtk`](https://bioinf.shenwei.me/csvtk) lets you manipulate tab-separated output of `gregex` at the command line; for example: 

```
$ python -m gregex 'Ab4GNb2(Ab4GNb4)Ma3' -o '|' -c | csvtk tab2csv -H -t | csvtk add-header -n Left,Match,Right | csvtk csv2md
Left           |Match    |Right
:--------------|:--------|:----
Ab4GNb2        |(Ab4GNb4)|Ma3
Ab4GNb2(Ab4GNb4|)        |Ma3


$ alias matches2md='csvtk tab2csv -H -t | csvtk add-header -n Left,Match,Right | csvtk csv2md'
$ python -m gregex 'Ab4GNb2(Ab4GNb4)Ma3' -o '|' -c | matches2md
Left           |Match    |Right
:--------------|:--------|:----
Ab4GNb2        |(Ab4GNb4)|Ma3
Ab4GNb2(Ab4GNb4|)        |Ma3


$ python -m gregex 'Ab4GNb2(Ab4GNb4)Ma3' -o '...' -c | csvtk tab2csv -H -t | csvtk add-header -n Left,Match,Right | csvtk csv2md
Left            |Match              |Right
:---------------|:------------------|:---------------
                |Ab4                |GNb2(Ab4GNb4)Ma3
                |Ab4GNb2            |(Ab4GNb4)Ma3
                |Ab4GNb2(Ab4GNb4)   |Ma3
                |Ab4GNb2(Ab4GNb4)Ma3|
Ab4             |GNb2               |(Ab4GNb4)Ma3
Ab4             |GNb2(Ab4GNb4)      |Ma3
Ab4             |GNb2(Ab4GNb4)Ma3   |
Ab4GNb2         |(Ab4GNb4)          |Ma3
Ab4GNb2         |(Ab4GNb4)Ma3       |
Ab4GNb2(        |Ab4                |GNb4)Ma3
Ab4GNb2(        |Ab4GNb4            |)Ma3
Ab4GNb2(Ab4     |GNb4               |)Ma3
Ab4GNb2(Ab4GNb4)|Ma3                |


$ python -m gregex 'Ab4GNb2(Ab4GNb4)Ma3' -o '_' -c | csvtk tab2csv -H -t | csvtk add-header -n Left,Match,Right | csvtk csv2md
Left            |Match              |Right
:---------------|:------------------|:---------------
                |Ab4                |GNb2(Ab4GNb4)Ma3
                |Ab4GNb2            |(Ab4GNb4)Ma3
                |Ab4GNb2(Ab4GNb4)   |Ma3
                |Ab4GNb2(Ab4GNb4)Ma3|
Ab4             |GNb2               |(Ab4GNb4)Ma3
Ab4             |GNb2(Ab4GNb4)      |Ma3
Ab4             |GNb2(Ab4GNb4)Ma3   |
Ab4GNb2         |(Ab4GNb4)          |Ma3
Ab4GNb2         |(Ab4GNb4)Ma3       |
Ab4GNb2(        |Ab4                |GNb4)Ma3
Ab4GNb2(        |Ab4GNb4            |)Ma3
Ab4GNb2(        |Ab4GNb4)           |Ma3
Ab4GNb2(        |Ab4GNb4)Ma3        |
Ab4GNb2(Ab4     |GNb4               |)Ma3
Ab4GNb2(Ab4     |GNb4)              |Ma3
Ab4GNb2(Ab4     |GNb4)Ma3           |
Ab4GNb2(Ab4GNb4 |)                  |Ma3
Ab4GNb2(Ab4GNb4 |)Ma3               |
Ab4GNb2(Ab4GNb4)|Ma3                |
```

## TODO

1. Migrate tests from the `dev` Jupyter notebook into `pytest` tests.
2. Add additional tests for code unique to `gregex.py` relative to the dev notebook (e.g. make sure parser recognizes every uncertainty-operator-free linear code expression you can find).
3. Create a clean demo notebook from the existing development notebook.
4. Setup `readthedocs` documentation.
5. Qualify imports in `gregex.py` to avoid polluting user namespace when `gregex` is imported as a module. 
6. Allow for distinct grammars to be loaded or swapped programmatically or specified via file (and supported through the CLI).
7. Add feature for stricter checking/enforcement of child ordering conventions.
8. Add support to the parser for uncertainty operators via a tool like `minikanren` or `z3`. Note that both directions will likely have limited support for Python 2.
9. Replace the parsing backend with something more efficient. There are many options here; ideally whatever is chosen should support more expressive grammars (e.g. left-recursive rules). 
10. Add pretty-printing support to s-expression conversion and make argument labels (=bond information) more explicit. For example, `NNa3(ANb4)Ab4GNb2(NNa6Ab4GNb4)Ma3(NNa3(ANb4)Ab4GNb3Ab4GNb2(NNa3(ANb4)Ab4GNb6)Ma6)Ma4GNb4(Fa6)GN`, when converted to an s-expression, should become something like one of these two examples below

```
(GN Fa6
    (GNb4 (Ma4 (Ma6 (GNb6 (Ab4 ANb4
                               NNa3))
                    (GNb2 (Ab4 (GNb3 (Ab4 ANb4
                                          NNa3)))))
               (Ma3 (GNb4 (Ab4 NNa6))
                    (GNb2 (Ab4 ANb4
                               NNa3))))))
```

```
(GN :a6 F
    :b4 (GN :a4 (M :a6 (M :b6 (GN :b4 (A :b4 AN
                                         :a3 NN))
                          :b4 (GN :b4 (A :b3 (GN :b4 (A :b4 AN
                                                        :a3 NN)))))
                   :a3 (M :b4 (GN :b4 (A :a6 NN))
                          :b2 (GN :b4 (A :b4 AN
                                         :a3 NN))))))
```

[`hy`](https://docs.hylang.org/en/stable/) plausibly has pretty-printing facilities that support this out-of-the-box.

