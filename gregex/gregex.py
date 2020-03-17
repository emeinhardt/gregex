from funcy import *
from itertools import product

from copy import deepcopy

from collections import OrderedDict
from json import dumps

import glypy
#from glypy.plot import plot
import glypy.io.linear_code

###########################################
# glypy imports and convenience functions #
###########################################

monosaccharides = glypy.monosaccharides

def gen(glycan):
    '''
    Abbreviation for
      `glypy.io.linear_code.to_linear_code(glycan)`
    '''
    return glypy.io.linear_code.to_linear_code(glycan)


def parse(linear_code_expression):
    '''
    Abbreviation for
      `glypy.io.linear_code.parse_linear_code(linear_code_expression)`
    '''
    return glypy.io.linear_code.parse_linear_code(linear_code_expression)


def parsePlot(linear_code_expression):
    '''
    Abbreviation for
      `plot(glypy.io.linear_code.parse_linear_code(linear_code_expression))`
    '''
    plot(glypy.io.linear_code.parse_linear_code(linear_code_expression))


#######################################
# defining some convenient categories #
#######################################


parentheses = ('(',')')

MS_codes = set()
has_no_code = set()

for each in list(monosaccharides):
    try:
        result = gen(monosaccharides.get(each))
        MS_codes.add(result)
        if result[-1] == '?':
            MS_codes.add(result[:-1])
    except Exception as e:
        has_no_code.add(each)

MS_codes.add('GN')
MS_codes.add('GNa')

SUs = MS_codes

bonds = set(map(str, range(1,10)))

bond_type_and_loc = set(map(partial(str_join, ''), 
                            product("ab?", 
                                    list(map(str, range(1,10))) + ['?'])))


SUs_with_bonds = set(map(partial(str_join, ''),
                         product(SUs,
                                 bonds)))


#######################################
# Bottom up parsing utility functions #
#######################################

def has_balanced_parens(s):
    '''
    Indicates whether the string in question has balanced parentheses.
    '''
    paren_stack = []
    for x in s:
        if x == '(':
            paren_stack.append(x)
        if x == ')':
            if len(paren_stack) == 0:
                return False
#             top = paren_stack[-1]
            paren_stack.pop()
    return len(paren_stack) == 0


def every_left_paren_has_a_right_paren(s):
    '''
    Indicates whether every open parenthesis has a matching close parenthesis.
    '''
    paren_stack = []
    for x in s:
        if x == '(':
            paren_stack.append(x)
        if x == ')':
#             if len(paren_stack) == 0:
#                 return False
#             top = paren_stack[-1]
            if len(paren_stack) > 0:
                paren_stack.pop()
    return len(paren_stack) == 0


def generate_subsequences(s, as_generator=False, with_contexts=False):
    '''
    Given a sequence of length n, generates all O(n^2) contiguous subsequences
    of s.

    If with_contexts is True, returns a sequence of 3-tuples:
        ...(left context, subsequence, right context)...
    '''
    c = with_contexts
    subseqs = (s[i:j] if not c else (s[:i], s[i:j], s[j:])
               for i in range(len(s))
               for j in range(i+1,len(s)+1))
    if as_generator:
        return subseqs
    else:
        return tuple(subseqs)


to_str = partial(str_join, '')
rev_to_str = compose(to_str, reversed, to_str)

def generate_prefixes(s, as_generator=False):
    l = len(s)
    prefs = (s[:i] for i in range(l+1))
    if as_generator:
        return prefs
    return list(prefs)


def generate_suffixes(s, as_generator=False):
    s_prime = to_str(reversed(s)) if type(s) == str else tuple(reversed(s))
    if as_generator:
        if type(s) == str:
            return (rev_to_str(each)
                    for each in generate_prefixes(s_prime, as_generator))
        else:
            return (tuple(reversed(each))
                    for each in generate_prefixes(s_prime, as_generator))
    else:
        if type(s) == str:
            return list(map(rev_to_str, generate_prefixes(s_prime, as_generator)))
        else:
            return list(map(lambda s: tuple(reversed(s)),
                            generate_prefixes(s_prime, as_generator)))


def deprefix(prefix, s, prefix_not_found_behavior='identity'):
    '''
    Left-inverse of sequence concatenation: returns p\s, where p = prefix.
    
    Behavior when p is not a prefix of s is governed by 
    `prefix_not_found_behavior`:
      - 'identity' -> s
      - 'empty'    -> '' if s is a string else the empty tuple
      - 'null'     -> None
    '''
    l_prefix = len(prefix)
    l_s = len(s)
    prefix_match = False
    quotient = ''
    
    if l_prefix <= l_s:
        prefix_match = s[:l_prefix] == prefix
    
    if prefix_match:
        quotient = s[l_prefix:]
        return quotient
    
    if prefix_not_found_behavior == 'identity':
        return s
    elif prefix_not_found_behavior == 'empty':
        return '' if type(s) == str else tuple()
    elif prefix_not_found_behavior == 'null':
        return None
    else:
        raise Exception("behavior must be 'identity', 'empty', or 'null': got {0}".format(prefix_not_found_behavior))


def desuffix(suffix, s, suffix_not_found_behavior='identity'):
    '''
    Right-inverse of sequence concatenation: returns s/f where f = suffix.

    Behavior when f is not a suffix of s is governed by
    `suffix_not_found_behavior`:
      - 'identity' -> s
      - 'empty'    -> '' if s is a string else the empty tuple
      - 'null'     -> None
    '''
    if type(s) == str:
        result = deprefix(rev_to_str(suffix),
                          rev_to_str(s),
                          suffix_not_found_behavior)
    if type(s) == str and result is not None:
        return rev_to_str(result)
    if type(s) != str:
        result = deprefix(tuple(reversed(suffix)),
                          tuple(reversed(s)),
                          suffix_not_found_behavior)
    if type(s) != str and result is not None:
        return tuple(reversed(result))
    return result


def circumfix(prefix, suffix, s):
    '''
    Returns
        prefix + s + suffix
    '''
    return prefix + s + suffix


def decircumfix(prefix, suffix, s, circumfix_not_found_behavior='identity'):
    '''
    Inverse of circumfix(prefix, suffix, u): returns u = (prefix\s)/suffix.

    Behavior when s does not have (prefix, suffix) as a circumfix is governed by
    `circumfix_not_found_behavior`:
      - 'identity' -> s
      - 'empty'    -> '' if s is a string else the empty tuple
      - 'null'     -> None
    '''
    has_circumfix = False

    deprefixed = deprefix(prefix, s, circumfix_not_found_behavior)

    if deprefixed is None:
        return None

    if circumfix_not_found_behavior == 'identity' and prefix != '':
        has_circumfix = deprefixed != s
    if not has_circumfix and circumfix_not_found_behavior == 'identity':
        return s

    decircumfixed = desuffix(suffix, deprefixed, circumfix_not_found_behavior)

    if circumfix_not_found_behavior == 'identity' and suffix != '':
        has_circumfix = decircumfixed != deprefixed
    if not has_circumfix and circumfix_not_found_behavior == 'identity':
        return s

    return decircumfixed


def infix(u, left_context, right_context):
    '''
    Returns 
        left_context + u + right_context
    '''
    return left_context + u + right_context


def deinfix(u, s, infix_not_found_behavior='identity'):
    '''
    Multi-inverse of infix(u, left, right): returns all (left, right)s s.t.
        s = left + u + right
    
    Behavior when s does not have u as an infix is governed by
    `infix_not_found_behavior`:
      - 'identity' -> s
      - 'empty'    -> '' if s is a string else the empty tuple
      - 'null'     -> None
    '''
    all_subsequences_with_contexts = generate_subsequences(s, 
                                                           as_generator=False, 
                                                           with_contexts=True)
    relevant_subsequences = tuple(filter(lambda lcr: lcr[1] == u,
                                         all_subsequences_with_contexts))
    just_contexts = tuple(map(lambda lcr: (lcr[0], lcr[2]),
                              relevant_subsequences))
    if len(just_contexts) == 0:
        if infix_not_found_behavior == 'identity':
            return s
        elif infix_not_found_behavior == 'empty':
            return ''
        elif infix_not_found_behavior == 'null':
            return None
    
    return just_contexts


def tokenizer(linear_code_expression):
    '''
    Given a linear code expression for a single molecule s, splits ('tokenizes')
    s into 'tokens' with each token consisting of one of
        - parentheses
        - a saccharide unit with bond information
        - a saccharide unit/monosaccharide
    where 'a monosaccharide with bond information' is either
        - a monosaccharide with just a bond type
        - a monosaccharide with a bond type and a bond location

    (Note that it is NOT the job of a tokenizer to enforce or check either
    syntactic or semantic restrictions.)
    '''
    categories = (parentheses, SUs_with_bonds, SUs)

    def right_greedy_tokenizer(lce, tokenized_result_so_far):
        '''
        Greedily (and recursively) tokenizes lce from the right.
        '''
        if lce == '':
            return tokenized_result_so_far
#         print('lce={0}'.format(lce))
        for suf in generate_suffixes(lce, True):
#             print('\tsuf={0}'.format(suf))
            for cat in categories:
#                 print('\t\tcat={0}'.format(cat))
                if suf in cat:
                    rest = desuffix(suf, lce)
#                     print('\t\trest={0}'.format(rest))
                    return right_greedy_tokenizer(rest, [suf] +
                                                         tokenized_result_so_far)
        e = 'Tokenized portion: {0}\nUntokenized remainder: {1}'
        raise Exception(e.format(tokenized_result_so_far, lce))
    return right_greedy_tokenizer(linear_code_expression, [])


########################################
# Conversion to JSON and s-expressions #
########################################

def destem(tokens):
    if tokens[-1] == ')':
        return tuple([''])
    stem = tuple(reversed(list(takewhile(lambda token: token not in '()',
                                         list(reversed(tokens))))))
    rest = tuple(reversed(list(dropwhile(lambda token: token not in '()',
                                         list(reversed(tokens))))))
    return stem, rest


def get_rightmost_branch(nonleftmost_subtrees_as_token_seq):
    s = nonleftmost_subtrees_as_token_seq
    end_of_first_branch = len(s)
    for i in range(len(s)-1,-1,-1):
        if s[i] == '(' and is_nonleftmost_branch(s[i:]):
            return s[i:]
    raise Exception('Does not have a well-formed rightmost branch:\n{0}'.format(nonleftmost_subtrees_as_token_seq))
    
    
def parse_nonleftmost_branches(nonleftmost_subtrees_as_token_seq):
    s = nonleftmost_subtrees_as_token_seq
    
    #fixme convert this tail-recursive function to an iterative function...
    def split_nonleftmost_branch_helper(remainder, acc):
        if len(remainder) == 0:
            return acc, remainder
        elif remainder[-1] != ')':
            return acc, remainder
        rightmost_branch = get_rightmost_branch(remainder)
        new_remainder = desuffix(rightmost_branch, remainder)
        new_acc = acc + [rightmost_branch]
        return split_nonleftmost_branch_helper(new_remainder, new_acc)
    
    grouped_branches, tail = split_nonleftmost_branch_helper(s, [])
    unwrapped_groups = tuple(map(unwrap_branch_tokens, grouped_branches))
    parsed_groups = tuple(map(parse_tokens, unwrapped_groups))
    parsed_tail = parse_tokens(tail)
    return parsed_groups + tuple([parsed_tail])


def unwrap_branch_tokens(branch_tokens):
    s = branch_tokens
    assert s[0] == '(' 
    assert s[-1] == ')'
    return s[1:-1]


def parse_subtrees(rest_tokens):
    parsed_nonleftmost_subtrees = parse_nonleftmost_branches(rest_tokens)
    result = parsed_nonleftmost_subtrees
    return result

    
def parse_exp(linear_code_expression, style='stem-and-subtrees'):
    '''
    Converts `linear_code_expression` formatted according to this BNF rule 
    (assuming leftward-ascending normal form):
      (exp) non-leftmost-branch* stem <- exp
    as a dictionary (with ordering within sequences going from left-to-right).
    '''
    tokens = tokenizer(linear_code_expression)
    return parse_tokens(tokens, style)
    

def parse_tokens(tokens, style='stem-and-subtrees'):
    stem, rest = destem(tokens)
    if rest == tuple():
        stem_and_subtrees = {'stem':tuple(reversed(stem)),
                             'subtrees':rest}
    else:
        try:
            subtrees = parse_subtrees(rest)
        except Exception as e:
            print("Tried to parse\n{0}\n as 'rest'".format(rest))
            raise e
        stem_and_subtrees = {'stem':tuple(reversed(stem)), 'subtrees':subtrees}
    if style == 'stem-and-subtrees':
        return stem_and_subtrees
    elif style == 'func-and-args':
        return stem_and_subtrees_to_func_and_args(stem_and_subtrees)
    elif style == 's-exp':
        return func_and_args_to_sexps(stem_and_subtrees_to_func_and_args(stem_and_subtrees), True)

    
def stem_to_func_and_args(stem):
    if len(stem) == 1:
        tree = OrderedDict()
        tree['func'] = stem[0]
        tree['args'] = tuple()
        return tree
#         return {'func': stem[0], 'args':tuple()}
    tree = OrderedDict()
    tree['func'] = stem[0]
    tree['args'] = tuple([stem_to_func_and_args(stem[1:])])
    return tree
#     return {'func':stem[0], 'args':tuple([stem_to_func_and_args(stem[1:])])}
    
    
def rightmost_leaf(fa_tree):
    t = fa_tree
    if len(t['args']) == 0:
        return t
    return rightmost_leaf(t['args'][-1])
    
    
def stem_and_subtrees_to_func_and_args(stem_and_subtrees):
    stem_fa = stem_to_func_and_args(stem_and_subtrees['stem'])
    my_rightmost_leaf = rightmost_leaf(stem_fa)
    my_rightmost_leaf['args'] = list(map(stem_and_subtrees_to_func_and_args,
                                         stem_and_subtrees['subtrees']))
    return stem_fa

def func_and_args_to_sexps(func_and_args, unwrap_leaves=False, separate_args=False):
    tree = func_and_args
    if len(tree['args']) == 0:
        return '({0})'.format(tree['func'])
    s = '({0}'.format(tree['func'])
    body = ' '.join(['{0}'.format(func_and_args_to_sexps(arg, unwrap_leaves, False))
                     for arg in tree['args']])
    s += ' ' + body + ')'
#     print('pre-unwrapping: {0}'.format(s))
    
    if unwrap_leaves:
        split_s = s.split(' ')
        unwrapped = map(lambda token: unwrap_branch_tokens(token )if token[0] == '(' and token[-1] == ')' else token,
                        split_s)
        joined_s = ' '.join(unwrapped)
        s = joined_s
#     print('unwrapped: {0}'.format(s))
    
    def extract_and_format_arglabels(token):
        is_leaf = token[0] != '('# and token[-1] != ')'
#         assert is_leaf, "{0} is not a leaf".format(token)
        if is_leaf:
            arg_label = token[-2:] if token[-1] != ')' else token[-3:-1]
            arg = token[:-2] if token[-1] != ')' else token[:-3]
            end_symbol = '' if token[-1] != ')' else ''
            assert arg_label in bond_type_and_loc, "{0} is not a recognized bond type\nOriginal token = {1}".format(arg_label, token)

            formatted_arg = ':' + arg_label + ' ' + arg + end_symbol
            return formatted_arg
        assert token[0] == '(', '{0} is neither a leaf nor a head label'.format(token)
        arg_label = token[-2:]
        arg = token[:-2]
        assert arg_label in bond_type_and_loc, "{0} is not a recognized bond type\nOriginal token = {1}".format(arg_label, token)

        formatted_arg = ':' + arg_label + ' ' + arg
        return formatted_arg
    
    if separate_args:
        #FIXME
        print('Currently unsupported.')
#         print('s = {0}'.format(s))
#         split_s = s.split(' ')
#         print('split_s = {0}'.format(split_s))
#         formatted = [split_s[0]] + list(map(extract_and_format_arglabels, split_s[1:]))
#         joined_s = ' '.join(formatted)
#         s = joined_s
#         print('separated s = {0}'.format(s))
        
    return s


##############################################
# Krambeck et al. 2009 ligand `...` operator #
##############################################

def is_ligand_match(linear_code_expression):
    '''
    Indicates whether `linear_code_expression` (in its entirety) matches (i.e.
    could be substituted with/for) `...`.
    '''
    s = linear_code_expression
    if len(s) == 0:
        return True
    return has_balanced_parens(s)


def get_ligand_matches(linear_code_expression, as_generator=False, with_contexts=False):
    '''
    Returns the nonempty substrings within `linear_code_expression` that match
    Krambeck et al's `ligand` uncertainty operator `...`.
    '''
    s = linear_code_expression
    subsequences = generate_subsequences(list(tokenizer(s)), True, with_contexts)
#     subsequences = generate_subsequences(s, True)
    if as_generator:
        return (is_ligand_match(subseq)
                for subseq in subsequences)
    else:
        return tuple(filter(is_ligand_match, subsequences))


####################################################
# Krambeck et al. 2009 continuation  `_`  operator #
####################################################

def is_continuation_match(linear_code_expression):
    '''
    Indicates whether `linear_code_expression` (in its entirety) matches (i.e.
    could be substituted with/for `_`).
    '''
    s = linear_code_expression
    if len(s) == 0:
        return True
    return every_left_paren_has_a_right_paren(s)


def get_continuation_matches(linear_code_expression, as_generator=False, with_contexts=False):
    '''
    Returns the nonempty substrings within `linear_code_expression` that match 
    Krambeck et al's `continuation` uncertainty operator `_`.
    '''
    s = linear_code_expression
    subsequences = generate_subsequences(list(tokenizer(s)), True, with_contexts)
#     subsequences = generate_subsequences(s, True)
    if as_generator:
        return (is_continuation_match(subseq) 
                for subseq in subsequences)
    else:
        return tuple(filter(is_continuation_match, subsequences))


###########################################################
# Krambeck et al. 2009 possible branch point `|` operator #
###########################################################


def is_possible_branch_point_match(linear_code_expression):
    '''
    Indicates whether `linear_code_expression` (in its entirety) matches (i.e.
    could be substituted with/for) either `` or `(...)`.
    '''
    s = linear_code_expression
    if len(s) == 0:
        return True
    left_edge_is_left_paren = s[0] == '('
    right_edge_is_right_paren = s[-1] == ')'
    center_matches_ligand = is_ligand_match(s[1:-1])
    return all([left_edge_is_left_paren, 
                right_edge_is_right_paren, 
                center_matches_ligand])


def get_possible_branch_point_matches(linear_code_expression, 
                                      as_generator=False, 
                                      with_contexts=False):
    '''
    Returns the nonempty substrings within `linear_code_expression` that match 
    Krambeck et al's `possible branch point` uncertainty operator `...`.
    '''
    s = linear_code_expression
    subsequences = generate_subsequences(list(tokenizer(s)), True, with_contexts)
#     subsequences = generate_subsequences(s, True)
    if as_generator:
        return (is_possible_branch_point_match(subseq) 
                for subseq in subsequences)
    else:
        return tuple(filter(is_possible_branch_point_match, subsequences))



