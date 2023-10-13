#!/usr/bin/env python

'''Filters HMMER search Output for best hit.

We ran HMMER search with all 6 amino acid reading frames translated
from nucleotide sequence. There is a low rate of more than one reading
frame per read finding a match and we want to keep only the best match.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: May 2023
License :: GNU GPLv3
Copyright 2023 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, random


def best_hits(query, bitscore, d, entry, dups):
    """ Filters the besthit based on bitscore """

    if query in d:
        dups += 1
        old_bitscore = float(d[query][0][1])
        if bitscore > old_bitscore:
            d[query] = [entry]

        elif bitscore == old_bitscore:
            d[query].append(entry)

    else:
        d[query] = [entry]

    return d, dups


def hmmer_search_filter(infile):

    d = {} # initialize dictionary for bitscore besthits
    dups = 0 # counter for number of duplicate matches for reads
    total = 0 # counter for total hmmer entries in file

    with open(infile, 'r') as f:
        # skip the hmm header
        for _ in range(15):
            next(f)
        # read throught the results
        for l in f:
            # remove whitespace and split fields
            X = l.lstrip().rstrip().split()
            # break the loop before the domain annotation output
            if X[1] == 'inclusion': break
            # sort out results
            total += 1 # add total
            # define fields
            query = '_'.join(X[8].split('_')[1:])
            bitscore = float(X[1])
            # track duplicate hits 
            d, dups = best_hits(query, bitscore, d, X, dups)

    print('Total number of entries in hmmer file:', total)
    print('Number of duplicate hmmer matches:', dups)

    return d


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--in_file',
        help='Please specify the HMMER search output file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--out_file',
        help='Please specify the output file name!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    # define input param
    infile = args['in_file']
    outfile = args['out_file']
    print('\n\nRunning Script...\n')
    filtered_best_hits = hmmer_search_filter(infile)

    # Write output file
    with open(outfile, 'w') as o:
        for k,v in filtered_best_hits.items():
            X = random.choice(v)
            # Fields are query, evaule, score, bias
            lineout = f'{k}\t{X[0]}\t{X[1]}\t{X[2]}\n'
            o.write(lineout)
        print(
            'Number of best hit entries written to new file:',
            len(filtered_best_hits), '\n\n'
            )

if __name__ == "__main__":
    main()
