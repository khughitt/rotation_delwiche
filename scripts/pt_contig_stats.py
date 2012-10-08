#!/usr/bin/env python2
#-*- coding:utf-8 -*-
"""
Trinity contig statistics
"""
import sys
import os
import re

def main():
    """Main"""
    for filepath in sys.argv[1:]:
        lengths = []
        
        # get contigs lengths
        for line in open(filepath):
            if not line.startswith(">"):
                continue
            length = int(re.search('len=([\d]*)', line).groups()[0])
            lengths.append(length)
        
        # print number of contigs and average contig length
        name = os.path.basename(filepath)
        mean = sum(lengths) / float(len(lengths))
        
        print("%s (n=%d mean=%f)" % (name, len(lengths), mean))

if __name__ == "__main__":
    sys.exit(main())