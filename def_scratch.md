# Vocabulary and conventions

**prefix, suffix** of a sequence *always* refer to a contiguous subsequence at the *left* or *right* edge of the sequence, respectively.
Note that a prefix or suffix need not be *proper*, and hence can include the entire sequence.

**precedes, follows** are always synonymous with *to the left of* or *to the right of*, respectively.

# Ligand operator

A ligand operator match must have balanced parentheses.

Therefore, for a substring of a glycan in linear code to match a ligand operator, the entire matching substring must be either
  - the empty string
  - any infix of a `stem`
  - a `non_main_branch`, 
    - ...possibly with a suffix of a `subexp` to its left 
      - = a suffix of a `substem` (= part of a chain), or 
      - = a complete `substem` (= complete chain) preceded by one or more complete(!) `non_main_branch`es themselves preceded by a (possibly empty) suffix of a `subexp`
    - ...possibly with zero or more complete(!) `non_main_branch`es followed by a (possibly empty) prefix of a `stem` to its right

# Continuation operator

A continuation operator match must have a matching right parenthesis for every left parenthesis, but not every right parenthesis needs to be matched. 

Therefore, for a substring of a glycan in linear code to match a continuation operator, the entire matching substring must either 
  - meet the conditions to match the ligand operator, *or*
  - *almost* match a `non_main_branch`, where the 'almost' is that the substring need not 
    - contain the left parenthesis. 
    - contain a complete `subexp`.
  - = the substring must contain a suffix of a `subexp` followed by `)`
    - ...possibly with zero or more complete(!) `non_main_branch`es followed by a (possibly empty) prefix of a `stem` to its right 
