import argparse
import os

my_desc = """Examine matches for uncertainty operators in a linear code expression.

Given
 - a linear code expression representing a single glycan
 - one of Krambeck et al. 2009's uncertainty operators ('...', '_', '|')
this returns a stream of lines (by default) representing information about
nonempty subsequences of the glycan that match the uncertainty operator.

If the context flag is not active AND no substitution arg is given, the stream
contains only the unique substrings of the glycan that match the operator.

If the context flag is active, each line of the stream contains three tab-
separated columns:
  left context\t match\t right context

If a substitution arg is provided, then this also returns an extra column for
each match (a boolean = 1 or 0) indicating whether both of the following 
conditions hold:
 - the substitution matches the operator
 - the expression resulting from the substitution is syntactically well-formed.
  - Currently, well-formedness just means that parentheses are balanced in
    the complete expression post-substitution.

If a filename is passed to the -x argument, then the script will write all data
to an xls-formatted file instead of stdout.
"""

parser = argparse.ArgumentParser(description=my_desc,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('lce', metavar ='LCE', 
                    type=str, nargs=1, 
                    help='a linear code expression containing no uncertainty operator tokens, or exactly one')
parser.add_argument('-o', '--operator', metavar='O',
                    type=str, nargs=1,
                    choices=(None, '...', '_', '|'),
                    help="one of Krambeck et al. 2009's three uncertainty operators, '...', '_', or '|'")
parser.add_argument('-c','--contexts',
                    action='store_true',
                    help='If active, then each match includes left and right context information.')
parser.add_argument('-s','--substitution', metavar='S',
                    type=str, nargs=1,
                    help='If provided, then the script checks whether the single token of a unique uncertainty operator in `LCE` can match `S` and whether the resulting linear code expression is well-formed.')
parser.add_argument('-x','--excel', metavar='X',
                    type=str, nargs=1,
                    help='If a filepath is provided via this arg, then output is written to an excel-formatted file at this location.')

args = parser.parse_args()
print(args)
