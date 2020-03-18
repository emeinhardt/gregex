import sys
import argparse

#print('foo')

#my_args = sys.argv
#for i, arg in enumerate(my_args):
#    print("{0}, {1}".format(i, arg))

my_desc = """Examine matches for uncertainty operators in a linear code expression.

Given
 - a linear code expression representing a single glycan
 - one of Krambeck et al. 2009's uncertainty operators ('...', '_', '|')
   this returns a stream of lines (by default) representing information about
   nonempty subsequences of the glycan that match the uncertainty operator.

If the context flag is not active AND no substitution arg is given, the stream
contains the unique substrings of the glycan that match the operator.

If the context flag is active, each line contains three tab-separated columns:
  left context	match	right context

If a substitution arg is provided, then this also returns an extra column for
each match (a boolean = 1 or 0) indicating whether both
 - the substitution matches the operator
 - the expression resulting from the substitution is syntactically well-formed.
  - Currently, well-formedness just means that parentheses are balanced in
    the complete expression post-substitution.

If a filename is passed to the -x argument, then the script will write all data
to an xls formatted file instead of stdout.
"""

parser = argparse.ArgumentParser(description=my_desc)
parser.add_argument('lce', metavar ='e', type=str, nargs=1, help='a linear code expression containing no uncertainty operator tokens, or exactly one')
parser.add_argument('--contexts', metavar='c', type=bool, nargs=1, help='If flag is active, then results ')

args = parser.parse_args()

