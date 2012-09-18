#!/usr/bin/env python2
#-*- coding:utf-8 -*-
"""
Generate E-Value histogram plots for each HMMER csv file specified

Usage:

    pt_evalue_histograms.py Mviride.csv Nmirabilis.csv
"""
import sys
import hmmer
import os

def main():
    """Main"""
    recarrays = [hmmer.parse_csv(x) for x in sys.argv[1:]]
    
    for i, filepath in enumerate(sys.argv[1:]):
        filename, ext = os.path.splitext(os.path.basename(filepath))
        title = "%s TAP domain" % filename
        hmmer.plot_evalue_histogram(recarrays[i], filename, title)

if __name__ == '__main__':
    sys.exit(main())
