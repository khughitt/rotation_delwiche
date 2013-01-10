TAP Classification Helper Scripts
=================================

This directory includes some short helpers scripts produced during the TAP
classification project.

Contents
--------
* **count_tap_domains.py** - Parses hmmsearch output and counts the number of 
TAP domains found for each species.
* **evalue_histograms.py** - Generates histograms of hmmsearch E-values.
* **hmmer.py** - Methods for working with HMMER3 output. Includes a method to
convert hmmsearch output to a [NumPy recarray](http://docs.scipy.org/doc/numpy/reference/generated/numpy.recarray.html)
as well as a method for generating a histogram of hmmsearch E-values.
* **pfam_comparison.py** - Compares output from a full Pfam-A.hmm search for a 
number of target speicies.
* **trinity_contig_stats.py** - Parses [Trinity](http://trinityrnaseq.sourceforge.net/)
FASTA output and prints number of contigs and mean contig length.


