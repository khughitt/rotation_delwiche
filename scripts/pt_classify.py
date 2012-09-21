#!/usr/bin/env python2
#-*- coding:utf-8 -*-
"""
Plant TAP domain classification
Keith Hughitt <khughitt@umd.edu>
2012/09/16

Classifies matched protein domains based on the rules described in 
Lang et al. (2010; doi: 10.1093/gbe/evq032)

@TODO 2012/09/20: deal with multiple classifications (prioritize TF 
classification)


@TODO: verify classification rules
@TODO: make sure no contigs are classified more than once
@TODO: very contigs that are ruled out

=> use E-value scores from HMMsearch output?

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
    for recarray in [hmmer.parse_csv(x) for x in sys.argv[1:]]:
        # list to store contigs in
        contigs = []
        
        # sort domain matches by contig id
        for contig_id in set(recarray['target_name']):
            domains = []
            
            # add unique domains associated with the contig id
            for row in recarray[recarray['target_name'] == contig_id]:
                if row['query_name'] not in [d.name for d in domains]:
                    domain = ProteinDomain(row['query_name'], row['Evalue'])
                    domains.append(domain)
            
            # create contig and add to the list
            contigs.append(Contig(contig_id, domains))

    # classify contigs
    for contig in contigs:
        for protein_family in protein_families:
            if protein_family.in_family(contig):
                print("%s: %s" % (contig.name, protein_family.name))
    
def init_classification_rules():
    "Initializes protein classification rules"""
    # @TODO: REPLACE legacy ids, add note...
    return [
        ProteinFamily("ABI3/VP1", "TF", ["B3"], ["AP2", "Auxin_resp", "WRKY"]),
        ProteinFamily("Alfin-like", "TF", ["Alfin-like"], ["Homeobox", "zf-TAZ", "PHD"]),
        ProteinFamily("AP2/EREBP", "TF", ["AP2"]),
        ProteinFamily("ARF", "TF", ["Auxin_resp"]),
        ProteinFamily("Argonaute", "TR", ["Piwi", "PAZ"]),
        ProteinFamily("ARID", "TF", ["ARID"]),
        ProteinFamily("AS2/LOB", "TF", ["DUF260"], ["bZIP_1", "bZIP_2", "HLH", "Homeobox"]),
        ProteinFamily("Aux/IAA", "TR", ["AUX_IAA"], ["Auxin_resp", "B3"]),
        ProteinFamily("BBR/BPC", "TF", ["DUF1004"]),
        ProteinFamily("BES1", "TF", ["DUF822"]),
        ProteinFamily("bHLH", "TF", ["HLH"]),
        ProteinFamily("bHSH", "TF", ["TF_AP-2"]),
        ProteinFamily("BSD domain containing", "PT", ["BSD"]),
        ProteinFamily("bZIP", "TF", [], ["HLH", "Homeobox"], alt_domains=["bZIP_1", "bZIP_2"]),
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
        ProteinFamily("Coactivator p15", "TR", ["PC4"]),
        ProteinFamily("CPP", "TF", ["CXC"]),
        ProteinFamily("CSD", "TF", ["CSD"]),
        ProteinFamily("CudA", "TF", ["STAT_bind", "SH2"]),
        ProteinFamily("DBP", "TF", ["DNC", "PP2C"]),
        ProteinFamily("DDT", "TR", ["DDT"], ["Homeobox", "Alfin-like"]),
        ProteinFamily("Dicer", "TR", ["DEAD", "Helicase_C", "Ribonuclease_3", "dsrm"], ["Piwi"]),
        ProteinFamily("DUF246 domain containing", "PT", ["DUF246"]),
        ProteinFamily("DUF296 domain containing", "PT", ["DUF296"]),
        ProteinFamily("DUF547 domain containing", "PT", ["DUF547"]),
        ProteinFamily("DUF632 domain containing", "PT", ["DUF632"]),
        ProteinFamily("DUF833 domain containing", "PT", ["DUF833"]),
        ProteinFamily("E2F/DP", "TF", ["E2F_TDP"]),
        ProteinFamily("EIL", "TF", ["EIN3"]),
        ProteinFamily("FHA", "TF", ["FHA"]),
        ProteinFamily("GARP_ARR-B", "TF", ["Response_reg"], ["CCT"], alt_domains=["G2-like_Domain", "Myb_DNA-binding"]),
        ProteinFamily("GARP_G2-like", "TF", ["G2-like_Domain"], ["Response_reg", "Myb_DNA-binding"]),
        ProteinFamily("GeBP", "TF", ["DUF573"]),
        ProteinFamily("GIF", "TR", ["SSXT"]),
        ProteinFamily("GNAT", "TR", ["Acetyltransf_1"], ["PHD"]),
        ProteinFamily("GRAS", "TF", ["GRAS"]),
        ProteinFamily("GRF", "TF", ["QLQ", "WRC"], []),
        ProteinFamily("HB", "TF", ["Homeobox"], ["EIN3", "KNOX1", "KNOX2", "HALZ", "bZIP_1"]),
        ProteinFamily("HB_KNOX", "TF", ["KNOX1", "KNOX2"]),
        ProteinFamily("HD-Zip", "TF", ["Homeobox"], alt_domains=["HALZ", "bZIP_1"]),
        ProteinFamily("HMG", "TR", ["HMG_box"], ["ARID", "YABBY"]),
        ProteinFamily("HRT", "TF", ["HRT"]),
        ProteinFamily("HSF", "TF", ["HSF_DNA-bind"]),
        ProteinFamily("IWS1", "TR", ["IWS1_C"], []),
        ProteinFamily("Jumonji", "TR", ["JmjC", "JmjN"], ["ARID", "GATA", "zf-C2H2", "Alfin-like"]),
        ProteinFamily("LFY", "TF", ["FLO_LFY"]),
        ProteinFamily("LIM", "TF", ["two_or_more_LIM"]),
        ProteinFamily("LUG", "TR", ["LUFS_Domain"]),
        ProteinFamily("MADS", "TF", ["SRF-TF"]),
        ProteinFamily("MBF1", "TR", ["MBF1"]),
        ProteinFamily("MED6", "TR", ["MED6"]),
        ProteinFamily("MED7", "TR", ["MED7"]),
        ProteinFamily("mTERF", "TF", ["mTERF"]),
        ProteinFamily("MYB", "TF", ["two_or_more_Myb_DNA-binding"], ["ARID", "G2-like_Domain", "Response_reg", "trihelix"]),
        ProteinFamily("MYB-related", "TF", ["Myb_DNA-binding"], ["ARID", "G2-like_Domain", "Response_reg", "trihelix", "two_or_more_Myb_DNA-binding"]),
        ProteinFamily("NAC", "TF", ["NAC_plant"]),
        ProteinFamily("NZZ", "TF", ["NOZZLE"]),
        ProteinFamily("OFP", "TR", ["DUF623"]),
        ProteinFamily("PcG_EZ", "TR", ["SANTA", "SET"]),
        ProteinFamily("PcG_FIE", "TR", ["FIE_clipped_for_HMM", "WD40"]),
        ProteinFamily("PcG_VEFS", "TR", ["VEFS-Box"], ["zf-C2H2"]),
        ProteinFamily("PHD", "TR", ["PHD"], ["Myb_DNA-binding", "Alfin-like", "ARID", "DDT", "Homeobox", "JmjC", "JmjN", "SWIB", "zf-TAZ", "zf-MIZ", "zf-CCCH"]),
        ProteinFamily("PLATZ", "TF", ["PLATZ"]),
        ProteinFamily("Pseudo ARR-B", "TF", ["CCT", "Response_reg"], ["tify"]),
        ProteinFamily("RB", "TF", ["RB_B"]),
        ProteinFamily("Rcd1-like", "TR", ["Rcd1"]),
        ProteinFamily("Rel", "TF", ["RHD"]),
        ProteinFamily("RF-X", "TF", ["RFX_DNA_binding"]),
        ProteinFamily("RRN3", "TR", ["RRN3"]),
        ProteinFamily("Runt", "TF", ["Runt"]),
        ProteinFamily("RWP-RK", "TF", ["RWP-RK"]),
        ProteinFamily("S1Fa-like", "TF", ["S1FA"]),
        ProteinFamily("SAP", "TF", ["STER_AP"]),
        ProteinFamily("SBP", "TF", ["SBP"]),
        ProteinFamily("SET", "TR", ["SET"], ["zf-C2H2", "CXC", "PHD", "Myb_DNA-binding"]),
        ProteinFamily("Sigma70-like", "TF", ["Sigma70_r2", "Sigma70_r3", "Sigma70_r4"]),
        ProteinFamily("Sin3", "TR", ["PAH"], ["WRKY"]),
        ProteinFamily("Sir2", "TF", ["SIR2"]),
        ProteinFamily("SOH1", "TR", ["SOH1"]),
        ProteinFamily("SRS", "TF", ["DUF702"]),
        ProteinFamily("SWI/SNF_BAF60b", "TR", ["SWIB"]),
        ProteinFamily("SWI/SNF_SNF2", "TR", ["SNF2_N"], ["AP2", "PHD", "zf_CCCH", "Myb_DNA-binding"]),
        ProteinFamily("SWI/SNF_SWI3", "TR", ["SWIRM"], ["Myb_DNA-binding"]),
        ProteinFamily("TAZ", "TF", ["zf-TAZ"]),
        ProteinFamily("TCP", "TF", ["TCP"]),
        ProteinFamily("TEA", "TF", ["TEA"]),
        ProteinFamily("TFb2", "TR", ["Tfb2"]),
        ProteinFamily("tify", "TF", ["tify"]),
        ProteinFamily("TRAF", "TR", ["BTB"], ["zf-TAZ"]),
        ProteinFamily("Trihelix", "TF", ["trihelix"]),
        ProteinFamily("TUB", "TF", ["Tub"]),
        ProteinFamily("ULT", "TF", ["ULT_Domain"]),
        ProteinFamily("VARL", "TF", ["VARL"]),
        ProteinFamily("VOZ", "TF", ["VOZ_Domain"]),
        ProteinFamily("Whirly", "TF", ["Whirly"]),
        ProteinFamily("WRKY", "TF", ["WRKY"]),
        ProteinFamily("zf_HD", "TF", ["ZF-HD_dimer"]),
        ProteinFamily("Zinc_finger, AN1 and A20 type", "TR", ["zf-AN1"], ["zf-C2H2"]),
        ProteinFamily("Zinc finger, MIZ type", "TF", ["zf-MIZ"], ["zf-C2H2"]),
        ProteinFamily("Zinc finger, ZPR1", "TR", ["zf-ZPR1"]),
        ProteinFamily("Zn_clus", "TF", ["Zn_clus"])
    ]
    
class Contig(object):
    """Class representing a single Contig, including all of it's matched
    domains.
    
    In cases where similar domains are encountered, the highest-scoring domain
    is kept, per the TAP classification guidelines.
    """
    def __init__(self, name, domains):
        self.name = name
        self.domains = domains
        
        # similar domains
        self.similar_domains = [
            set(["NF-YB", "NF-Y3", "CCAAT-Dr1_Domain"]),
            set(["PHD", "Alfin-like"]),
            set(["G2-like_Domain", "Myb_DNA-binding"]),
            set(["GATA", "zf-Dof"])
        ]
        
        # filter out similar domains and create a set with just the domain
        # names for comparison
        self._filter_similar()         

    def get_domain_names(self):
        """Returns a list of the matching domain names to use during 
        classification"""
        return set([d.name for d in self.domains])
    
    def _filter_similar(self):
        """Check for similar domain matches and keep only the closest hit"""
        for similar in self.similar_domains:
            domain_names = self.get_domain_names()
            
            intersection = domain_names.intersection(similar)
            
            # if more than one similar domain exists
            if len(intersection) > 1:
                # find top hit
                top_hit = sorted(self.domains, 
                                 key=lambda domain: domain.evalue,
                                 reverse=True).pop()

                # remove all other domains
                lower_hits = intersection.difference([top_hit.name])
                
                self.domains = filter(lambda d: d.name not in lower_hits, 
                                      self.domains)
            
class ProteinDomain(object):
    """Class representing a single protein domain"""
    def __init__(self, name, evalue):
        """Creates a new ProteinDomain instance"""
        self.name = name
        self.evalue = evalue

class ProteinFamily(object):
    """Class representing a Protein Family classfication rule"""
    def __init__(self, name, type_, requires, forbids=None, alt_domains=None):
        """Creates a new ProteinFamily instance"""
        self.name = name
        self.type = type_
        self.requires = set(requires)
        
        # Forbidden domains
        if forbids is not None:
            self.forbids = set(forbids)
        else:
            self.forbids = set()
            
        # Alternate domains (one of two must be present)
        if alt_domains is not None:
            self.alt_domains = set(alt_domains)
        else:
            self.alt_domains = set()
    
    def in_family(self, contig):
        """Checks to see whether a contig belongs to the protein family.
        
        Parameters:
        -----------
        contig : set
            Set of protein domains associated with the contig
        """
        contig_domains = contig.get_domain_names()
        
        return (self.requires.issubset(contig_domains) and 
                self.forbids.isdisjoint(contig_domains) and
                (len(self.alt_domains) == 0 or 
                 len(self.alt_domains.intersection(contig_domains)) > 0))
    
if __name__ == '__main__':
    #sys.exit(main())
    rules = main() # TEMP: Debugging


