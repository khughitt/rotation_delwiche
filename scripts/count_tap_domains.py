#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Parses HMMsearch table output and generate a CSV file of TAP domain counts for
each species.

Keith Hughitt <khughitt@umd.edu>
2012/09/15

Usage:
    count_tap_domains.py ../input/hmmsearch_taps/*
"""
import sys
import os
import csv
import itertools
import numpy as np
import hmmer

def main():
    """Main application"""
    # parse csv files for each species
    recarrays = [hmmer.parse_csv(x) for x in sys.argv[1:]]

    # protein domain counts csv
    write_counts_csv(recarrays, '../csv/domains/tap_domain_frequencies.csv')

def write_counts_csv(recarrays, output_file):
    """Creates a new csv file of the domain counts for each species."""
    # write results as CSV
    writer = csv.writer(open(output_file, 'wt'))

    # write header row
    header = ["domain"]
    for filepath in sys.argv[1:]:
        header.append(os.path.basename(filepath))
    writer.writerow(header)

    # get a list of all of the domains found
    for domain in get_domains(recarrays):
        row = [domain]
        for species in recarrays:
            row.append(count_matches(species, 'query_name', domain))
        writer.writerow(row)

def count_matches(recarray, field, value):
    """Returns the number of occurances of a given value for the specified
    field"""
    return np.count_nonzero(recarray[field] == value)

def get_domains(recarrays):
    """Returns a unique list of matched protein domains"""
    # create a list of the domain columns for each recarray
    matched_domains = [x['query_name'] for x in recarrays]

    # flatten into a single list
    flattened = list(itertools.chain(*matched_domains))

    return np.unique(flattened)

if __name__ == "__main__":
    sys.exit(main())
