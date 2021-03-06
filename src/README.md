# Code Structure

1. `scan.cpp`: The goal of this step is to build the prefix free parsing.
It first parses the reference, then extends the read with parts of the reference so that they start and end with a triggering sub-string (if possible), and then parses this extended read. To avoid going through the input file twice, the parse is first stored with a random hash for each phrases in a intermediate file `file.parse_old`. Before the start of a read there is a separator (which is PRIME+1) followed by the position where the extended read starts in the parse.
So that later we can avoid indexing characters belonging to the extensions and not to the original input (genome or read), we store for each phrase in the parse the limits of what needs to be index. This is stored in `file.limits`.
After the input (reference and reads) has been entirely parsed, the dictionary is sorted in **reverse lexicographic order** and outputted in the file `file.dict` and the number of occurrences of each word is stored in `file.occ`.
Finally the parse in `file.parse_old` is mapped to match the sorted dictionary. We replace the hash of a phrase by it's rank in reverse lexicographic order. If a word is the first word of the dictionary, it's hash is mapped to 1. The separator is now number_of_words_in_the_dictionary+1 and is put at the beginning of the file where we store this new and final parse : `file.parse`.
    - Input: `file`
    - Output: `file.parse_old`, `file.dict`, `file.occ`, `file.parse`,`file.limits`
2. `bwtparse.cpp`: The goal of this step is to build the inverted list of the parse (with regard to the alphabet of the parse), to be used for the final step of the construction.
From the `file.parse` we obtain the tree structure of the reads attached to the reference, then we apply a doubling algorithm to trace back to the root of the tree (the start of reference). This way each letter is sorted based on it's context (in the read followed by the reference) with regard to the **reverse lexicographic order**. Then we modify our structure of the tree to get list of the children of a node to be able to construct the BWT and to get the list of word that follow a given word, we save this to `file.full_children`.So we can access it easily in the next phase we permute the limits stored in `file.limits` according to the XBWT of the parse, this is stored in `file.xbwt_limits`.  Finally we construct the inverted list that we output to `file.ilist`.
    - Input: `file.parse`, `file.occ`, `file.limits`
    - Output: `file.ilist`, `file.full_children`, `file.xbwt_limits`
3. `pfbwt.cpp`: This last steps builds the BWT of the input based on the dictionary (more precisely the suffix array of all inverse words in the dictionary), the inverted list of the parse and the list that contains, for a given word, the list of words that follow it.
First the suffix array of all the inverse words (the fact that they are inverse allow us to use a suffix array instead of a prefix array) in the dictionary is computed by using `gSACA-K`. Then taking each suffix in order of the sa, we add the corresponding characters according to three case:
    - If the length of the suffix is less that `w` then this character is accounted for in an other word. We can move on to the next suffix.
    - If the suffix is actually a full word, then the characters we need to output are actually the ones at the end (position `w+1` from the right) of the reversed words that follow this word. We obtain them by the dictionary in `file.dict` and `file.full_children`.
    - In the last case, this means we are inside a word and need to be carefull because another word might have the same suffix. We use the longest comon prefix computed at the same time as the suffix array on the dictionary to obtain all the words that have this suffix and ouput the rights characters by two techniques. First, if all those words have the same character that precedes the suffix, then we can just ouput the character according to the number of occurences (in `file.occ`) of each of those words. Else, this is the hardest case and we have to know which character to output when according to their context given by the inverted list in `file.ilist`.

  For all those cases before adding the character we check if it is part of the limits.
    - Input: `file.dict`, `file.occ`, `file.ilist`, `file.full_children`, `file.xbwt_limits`
    - Output: `file.bwt`

## Description of the intermediate files

- `file.parse_old`
    - contains a first a parse of the genome, separator (PRIME+1) and a sequence of triples, starting position of the parse of the parse of the extended read as 64 bit unsigned int in the parse of the reference, parse of the extended read as 64 bit ints, and a separator (PRIME+1) as 64 bit int. Positions are 0-based, the only difference wit `file.parse` is how we refer to the characters, here they are hashes of the word and in `file.parse` they correspond to the rank in colexicographic order. This diference also implies thediference of choice of separator.
- `file.parse`
  - containing the parse P of the reference and the reads with each word identified with its 1-based lexicographic rank (ie its position in D). The files starts by the separator, which is the number of words in the dictionionary (d) + 1, after that we have the parse of the reference and then for every read, a separator, the starting position of the parse of the read in the parse of the reference
 We assume the number of distinct words is at most 2^32-1.
  - Note that the words in the parsing overlap by wsize
  - 1-based for the letters of the parse
  - 0-based for the starting position of the parse of the extended reads in the parse of the reference.
  - size: 4p bytes (32 bit ints)
- `file.dict`
  - containing the reversed dictionary words sorted in colexicographic order (the non reversed word) with a 0x1 at the end of each word and a 0x0 at the end of the file.
  - the smallest dictionary word is d0 ="" that occurs once (but not is not encoded in the dictionary file)
  - first (imaginary) word is a dummy word therefore the parse identifies words by 1-based lexicographic rank
  - size: |D| + d + 1 (|D| is the sum of the word lengths)
- `file.occ`
  - the number of occurrences of each word in lexicographic order. We assume the number of occurrences of each word is at most 2^32-1
  - size: 4d bytes (32 bit ints)
- `file.limits`
  - for each phrase in the parse, the first and last character to index in that phrase.
  - 8p bytes
- `file.extended_input`
  - This is an optional file that can be obtained with `scan_extended.x` which outputs the extended reads and there starting position in the reference. The extended read do not contain the first `w` characters (the trigering substring) that overlaps with the reference. This is mostly for debugging purposes.
- `file.ilist` (inverted list occurrence)
  - For each dictionary word in lexicographic order, the list of BWT positions where that word appears (ie i \in ilist(w) <=> BWT[i]=w). There is also an entry for the empty word which is not in the dictionary but is assumed to be the smallest word.
  - size: 4p+4
  - `file.xbwt_limits`
    - for each phrase in the xbwt of parse, the first and last character to index in that phrase.
    - 8p bytes
- `file.full_children` This contains for every word, the list of word that follows it in the entire tree structure, including for the empty word. The limit between each word is marked by a separator (still d + 1) and this separator is the first character of the file.
  - size: 4d bytes (32 bit ints)
- `file.bwt` The resulting XBWT (containing exactly the character in the original input, the read extensions are not indexed)
