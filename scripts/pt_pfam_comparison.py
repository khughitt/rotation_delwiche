#!/usr/bin/env python2
#-*- coding:utf-8 -*-
"""Prints some basic statistics about full Pfam-A.hmm search results for 
a number of target species."""
import sys
import os
import csv
import hmmer

def main():
    """Main"""
    targets = {}
    
    # read in files
    for filepath in sys.argv[1:]:
        # species name
        name = os.path.splitext(os.path.basename(filepath))[0]
        
        # parse HMMER output and store in data dict
        recarray = hmmer.parse_csv(filepath)
        targets[name] = set(recarray['query_name'])
        
        # output basic statistics
        print("(%s) # Domains: %d (%d unique)" % (name, recarray.shape[0],
                                                  len(targets[name])))
        
    # pairwise comparisons
    write_correlation_csv(targets)
    write_correlation_csv(targets, True)
    

        
def write_correlation_csv(targets, normalize=False):
    """write csv"""
    # choose a file to write to
    if normalize:
        suffix = "_normalized"
    else:
        suffix = ""

    writer = csv.writer(open("../csv/domains/pfam_domains%s.csv" % (suffix), 'wt'))
    writer.writerow([None] + targets.keys())
    
    i = 0
    for species1, set1 in targets.items():
        row = [species1]
        j = 0
        for species2, set2 in targets.items():
            # normalized matrices are symmetric
            if normalize and j > i:
                row.append(None)
                continue
            
            # calculate union and intersections
            intersection = float(len(set1.intersection(set2)))
            union = float(len(set1.union(set2)))
            
            # if at least some domains are shared between species, calculate
            # a correlation
            numerator = intersection
            
            if normalize:
                denominator = union
            else:
                denominator = len(set1)
                
            # check for non-zero denominator and compute correlation
            if denominator > 0:
                correlation = numerator / denominator
            else:
                correlation = 0                    
                
            row.append(correlation)
            j += 1

        writer.writerow(row)
        i += 1

if __name__ == "__main__":
    result = main()