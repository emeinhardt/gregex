from funcy import *
import gregex

import argparse
#import os
from collections import OrderedDict
import csv

my_desc = """Examine matches for uncertainty operators in a linear code expression.

Given
 - a linear code expression representing a single glycan
 - one of Krambeck et al. 2009's uncertainty operators ('...', '_', '|')
this returns a stream of lines (by default) representing information about
nonempty subsequences of the glycan that match the uncertainty operator.

If the context flag (-c) is not active AND no substitution arg is given, the 
stream contains only the unique substrings of the glycan that match the 
operator.

If the context flag is active, each line of the stream contains three tab-
separated columns:
  left context\t match\t right context

If a substitution arg is provided (-s), then this also returns an extra column 
for each match (a boolean = 1 or 0) indicating whether both of the following 
conditions hold:
 - the substitution matches the operator
 - the expression resulting from the substitution is syntactically well-formed.
  - Currently, well-formedness just means that parentheses are balanced in
    the complete expression post-substitution.

If the column name flag (-n) is active, then the output will incldue a column 
header line before data. 

If a filename is passed to the -x argument, then the script will write all data
to an xls-formatted file instead of stdout.

If the verbose flag (-v) is active, information will be printed to stdout 
about calculation.
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
parser.add_argument('-n','--namecolumns',
                    action='store_true',
                    help='If active, then output will include a column header line')
parser.add_argument('-x','--excel', metavar='X',
                    type=str, nargs=1,
                    help='If an operator is provided via -o AND a filepath is provided via this arg, then output is written to an excel-formatted file at this location.')
parser.add_argument('-v','--verbose',
                    action='store_true',
                    help='If active, then prints extra information to stdout')

args = parser.parse_args()
if args.verbose:
    print(args)

with_context = args.contexts
to_excel_fp = args.excel[0] if args.excel is not None else None
lce = args.lce[0]
operator = args.operator
substitution = args.substitution
verbose = args.verbose
colnames = args.namecolumns

if substitution is not None and len(substitution) > 0:
    sub = substitution[0]
else:
    sub = None
if operator is not None and len(operator) > 0:
    op = operator[0]
else:
    op = None

bool_mapper = {True:1, False:0}

if sub is not None and op is None:
    if verbose:
        print('Substitution and linear code expression provided.\nChecking if substitution is valid...')
    if colnames:
        cols = ('valid_sub?',)
        col_string = str_join('\t', cols)
        print(col_string)
    print(gregex.check_match(lce, sub, verbose))
else:
    results = gregex.analyze_matches(lce, op, sub, with_context, verbose)
    if colnames:
        if with_context:
            cols = ['left_context','match', 'right_context']
        else:
            cols = ['match']
        if sub is not None:
            cols += ['valid_sub?']
        cols = tuple(cols)
        col_string = str_join('\t', cols)
        
    if to_excel_fp is None:
        if colnames:
            print(col_string)
        columnify = lambda match_result: str_join('\t', match_result)
        for result in results:
            if type(result) != tuple:
                print(result)
            else:
                #print(result)
                print(columnify(result))
    else:
        results_dict = list(map(OrderedDict,
                                map(lambda m: zip(cols, m),
                                    results)))
        with open(to_excel_fp, 'wb') as csv_file:
            excel_writer = csv.DictWriter(csv_file, fieldnames=cols)
            if colnames:
                excel_writer.writeheader()
            excel_writer.writerows(results_dict)

