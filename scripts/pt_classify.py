#!/usr/bin/env python2
#-*- coding:utf-8 -*-
"""
Plant TAP domain classification
Keith Hughitt <khughitt@umd.edu>
2012/09/16

Classifies matched protein domains based on the rules described in 
Lang et al. (2010; doi: 10.1093/gbe/evq032)

Usage:
------
pt_classify.py hmmertbl.csv hmmertbl2.csv...

References:
-----------
 http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2997552/
"""
import sys
import hmmer

def main():
    """Main application body""" 
    # initialize ProteinFamily instances for each set of classification rules
    protein_families = init_classification_rules()
    
    # read in hmmer tables
    for species in [hmmer.parse_csv(x) for x in sys.argv[1:]]:
        # create list of domains for each contig
        pass
        
    # TEMP: Debugging
    return protein_families
    
def init_classification_rules():
    "Initializes protein classification rules"""
    # @TODO: REPLACE legacy ids, add note...
    return [
        ProteinFamily("ABI3/VP1", "TF", ["B3"], ["AP2", "Auxin_resp", "WRKY"]),
        ProteinFamily("Alfin-like", "TF", ["Alfin-like"], ["Homeobox", "zf-TAZ", "PHD"]),
        ProteinFamily("AP2/EREBP", "TF", ["AP2"]),
        ProteinFamily("ARF", "TF", ["Auxin_resp"]),
        ProteinFamily("Argonaute", "TF", ["Piwi", "PAZ"]),
        ProteinFamily("ARID", "TF", ["ARID"]),
        ProteinFamily("AS2/LOB", "TF", ["DUF260"], ["bZIP_1", "bZIP_2", "HLH", "Homeobox"]),
        ProteinFamily("Aux/IAA", "TF", ["AUX_IAA"], ["Auxin_resp", "B3"]),
        ProteinFamily("BBR/BPC", "TF", ["DUF1004"]),
        ProteinFamily("BES1", "TF", ["DUF822"]),
        ProteinFamily("bHLH", "TF", ["HLH"]),
        ProteinFamily("bHSH", "TF", ["TF_AP-2"]),
        ProteinFamily("BSD domain containing", "TF", ["BSD"]),
        ProteinFamily("bZIP1", "TF", ["bZIP_1"], ["HLH", "Homeobox"]),
        ProteinFamily("bZIP2", "TF", ["bZIP_2"], ["HLH", "Homeobox"]),
        ProteinFamily("C2C2_CO-like", "TF", ["CCT", "zf-B_box"], ["GATA", "tify", "PLATZ"]),
        ProteinFamily("C2C2_Dof", "TF", ["zf-Dof"], ["GATA"]),
        ProteinFamily("C2C2_GATA", "TF", ["GATA"], ["tify", "zf-Dof"]),
        ProteinFamily("C2C2_YABBY", "TF", ["YABBY"]),
        ProteinFamily("C2H2", "TF", ["zf-C2H2"], ["zf-MIZ"]),
        ProteinFamily("C3H", "TF", ["zf-CCCH"], ["AP2", "SRF-TF", "two_or_more_Myb_DNA-binding", "zf-C2H2"]),
        ProteinFamily("CAMTA", "TF", ["CG-1", "IQ"]),
        ProteinFamily("CCAAT_Dr1", "TF", ["CCAAT-Dr1_Domain"], ["NF-YB", "NF-YC"]),
        ProteinFamily("CCAAT_HAP2", "TF", ["CBFB_NFYA"], ["bZIP_1", "b_ZIP2"]),
        ProteinFamily("CCAAT_HAP3", "TF", ["NF-YB"], ["CCAAT-Dr1_Domain", "NF-YC"]),
        ProteinFamily("CCAAT_HAP5", "TF", ["NF-YC"], ["CCAAT-Dr1_Domain", "NF-YB"]),
        ProteinFamily("Coactivator p15", "TF", ["PC4"]),
        ProteinFamily("CPP", "TF", ["CXC"]),
        ProteinFamily("CSD", "TF", ["CSD"]),
        ProteinFamily("CudA", "TF", ["STAT_bind", "SH2"]),
        ProteinFamily("DBP", "TF", ["DNC", "PP2C"]),
        ProteinFamily("DDT", "TF", ["DDT"], ["Homeobox", "Alfin-like"]),
        ProteinFamily("", "TF", [], []),
        ProteinFamily("", "TF", [], []),
        ProteinFamily("", "TF", [], []),
        ProteinFamily("", "TF", [], []),
        ProteinFamily("", "TF", [], []),
        ProteinFamily("", "TF", [], []),
        ProteinFamily("", "TF", [], []),
        ProteinFamily("", "TF", [], []),
        ProteinFamily("", "TF", [], []),
        ProteinFamily("", "TF", [], []),
        ProteinFamily("", "TF", [], []),
        ProteinFamily("", "TF", [], []),
        ProteinFamily("", "TF", [], []),
        ProteinFamily("", "TF", [], [])
    ]

class ProteinFamily(object):
    """Class representing a Protein Family classfication rule"""
    def __init__(self, name, type_, requires, forbids=None):
        """Creates a new ProteinFamily instance"""
        self.name = name
        self.type = type_
        self.requires = set(requires)
        
        if forbids is not None:
            self.forbids = set(forbids)
        else:
            self.forbids = set()
    
    def in_family(self, contig):
        """Checks to see whether a contig belongs to the protein family.
        
        Parameters:
        -----------
        contig : set
            Set of protein domains associated with the contig
        """
        return (self.requires.issubset(contig) and self.forbids.isdisjoint(contig))
    
if __name__ == '__main__':
    #sys.exit(main())
    rules = main() # TEMP: Debugging


