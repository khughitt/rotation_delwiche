#!/usr/bin/env python2
#-*- coding:utf-8 -*-
"""
HMMER2GO

keith hughitt <khughitt@umd.edu>
2012/10/07

Parses HMMER3 output tables and maps matched Pfam domains to GO terms.

Usage:
    hmmer2go.py gene_ontology.1_2.obo hmmer_results.csv
"""
import os
import re
import sys
import csv
import hmmer
from goatools import obo_parser

def analyze_hmmer_table(go_terms, pfam2go, hmmer_table, go_level=1):
    """Analyzes a given HMMER3 search table for Pfam/GO term matches.
    
    Arguments
    ---------
    go_terms : goatools.obo_parser.GODag
        A dictionary of GO terms
    pfam2go : dict
        Pfam/GO mapping
    hmmer_table : str
        filepath to a HMMER3 search output table
    go_level : int
        (Optional) the GO level to summarize
    """
    summary = {
        "unknown": 0
    }

    # parse HMMER output and store in data dict
    contigs = hmmer.parse_csv(hmmer_table)
    
    # iterate through HMMER output
    for contig in contigs:
        pfam_domain = contig[2]
        
        # check to see if domain is in Pfam2go
        if pfam_domain not in pfam2go:
            summary['unknown'] += 1
            continue
        
        # otherwise get a list of the associated GO terms
        terms = pfam2go[pfam_domain]
        
        # iterate through GO terms for the contig
        for t in terms:
            go_term = go_terms[t[0]]
    
            # find GO category at the desired level
            node = go_term
            
            for i in range(go_term.level - go_level):
                node = node.parents[0]
                
            # add to summary table
            category = node.name
            
            if not category in summary:
                summary[category] = 0
                
            summary[category] += 1
    
    return summary

def load_pfam2go(source=None):
    """Map HMMER3 results to Gene Ontology (GO) terms"""
    # download pfam2go mapping if none specified
    if source is None:
        import urllib2
        url = "http://www.geneontology.org/external2go/pfam2go"
        raw = urllib2.urlopen(url).read().split("\n")
    else:
        # otherwise read in local file
        fp = open(source, 'r')
        raw = fp.readlines()
        fp.close()
        
    # parse and store in a python dict
    mapping = {}
    
    # pfam2go entry regex
    regex = r"Pfam:(\w+) ([\w\-]+) > (GO:[\w\- \(\)\[\]\/\+>,':]+) ; (GO:\d+)"
    
    for line in raw:
        # skip comments
        if line.startswith("!") or line == '':
            continue
        
        # parse a single entry in pfam2go
        # ex: Pfam:PF00060 Lig_chan > GO:membrane ; GO:0016020
        match = re.match(regex, line)
        pfam_id, pfam_name, go_name, go_id = match.groups() 
        
        if pfam_name not in mapping:
            mapping[pfam_name] = []
            
        mapping[pfam_name].append((go_id, go_name))
        
    return mapping

def write_csv(results, filepath, norm=None):
    """Write GO analysis results to CSV"""
    # get a list of all GO categories used
    categories = set([])
    for data in results.values():
        categories = categories.union(data.keys())

    # if no normalization requested do not scale results
    if norm is None:
        scale_factors = {k: 1 for k in results.keys()}
    else:
        # otherwise divide norm value by the number of Pfam domains used
        scale_factors = {}
        
        for species, data in results.items():
            # get total number of domains
            total = 0
            
            for term, value in data.items():
                total += value
                
            # scale factor = NORM / NUM DOMAINS
            scale_factors[species] = norm / float(total)
        
    # write results to CSV
    writer = csv.writer(open(filepath, 'wt'))
    writer.writerow(["category"] + results.keys())
    
    for category in categories:
        row = [category]
        for species, data in results.items():
            # add category count if it exists for the species
            if category in data:
                row.append(int(round(data[category] * scale_factors[species])))
            else:
                # otherwise add 0
                row.append(0)
        # add row to CSV
        writer.writerow(row)
    

def parse_hmmer_table(filepath):
    """Parses a hmmer table and outputs a summary of the Pfam GO associations"""
    
if __name__ == "__main__":
    go_input = "../input/gene_ontology.1_2.obo"
    pfam2go_file = "../input/pfam2go.txt"
    
    # load Gene ontology
    print("Loading Gene Ontology...")
    go_terms = obo_parser.GODag(go_input)
    #go_terms = obo_parser.GODag(sys.argv[1])

    # load pfam2go
    print("Loading Pfam2go...")
    pfam2go = load_pfam2go(pfam2go_file)
    
    # map GO terms
    results = {}
    
    for filepath in sys.argv[1:]:
        name = os.path.splitext(os.path.basename(filepath))[0]
        results[name] = analyze_hmmer_table(go_terms, pfam2go, filepath)
    
    # output results to CSV
    write_csv(results, "../csv/go_analysis.csv")
    write_csv(results, "../csv/go_analysis_normed.csv", norm=5e4)

