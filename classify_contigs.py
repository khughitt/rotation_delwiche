#!/usr/bin/env python2
#-*- coding:utf-8 -*-
"""
Plant TAP domain classification
Keith Hughitt <khughitt@umd.edu>
2012/09/21

Classifies matched protein domains based on the rules described in 
Lang et al. (2010; doi: 10.1093/gbe/evq032)

@TODO: verify classification rules
@TODO 2012/09/23: verify redundant classification removal step; sort 
classifications by family
@TODO 2012/09/24: generate summary statistics for domain matches
@TODO 2012/09/24: summarize overlap in families between species (pairwise 
similarity?)

Usage:
------
classify_contigs.py hmmertbl.csv hmmertbl2.csv...

References:
-----------
 http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2997552/
"""
import sys
import csv
import os
import numpy as np
from scripts import hmmer

def main():
    """Main application body""" 
    # initialize ProteinFamily instances for each set of classification rules
    protein_families = init_classification_rules()
    
    # create a master dictionaries to keep track of results
    results = {}
    domains = {}
    
    # create necessary directories if they do not already exist
    create_output_dirs()
    
    # read in hmmer tables
    for filepath in sys.argv[1:]:
        classify_species(filepath, results, domains, protein_families)
        
    # generate a dictionary representing the diversity (unique) tap complements
    # for each target species
    tap_diversity = get_tap_diversity(results)
        
    # output summary files
    write_classification_summary(results, tap_diversity)
    
    #
    # TODO: Pickle results and break up into smaller scripts.
    #    
    write_domain_summary_csv(domains)
    write_tap_correlation_csv(tap_diversity)
    write_tap_correlation_csv(tap_diversity, normalize=True)
    write_extended_tap_correlation_csv(tap_diversity)
    
    # output phylip character matrices
    write_phylip_char_matrix(tap_diversity, protein_families, "cyame")
    write_phylip_char_matrix(tap_diversity, protein_families, "cyame", ["TF"], 
                             "infile_TF")
    write_phylip_char_matrix(tap_diversity, protein_families, "cyame", ["TR"], 
                             "infile_TR")

    return results

def classify_species(filepath, results, domains, protein_families):
    """Classifies a species using the results from hmmsearch"""
    recarray = hmmer.parse_csv(filepath)
    
    # output base filename
    target = os.path.splitext(os.path.basename(filepath))[0]
           
    # Load contigs with matching domains
    contigs = load_contigs(target, recarray, domains)

    # open csv file for output
    filename =  "classification_%s.csv" % target
    fp = open(os.path.join("../csv/classifications", filename), 'wt')
    csv_writer = csv.writer(fp)
    csv_writer.writerow(['contig', 'family', 'type'])

    # classify contigs
    classifications = classify_contigs(contigs, protein_families)

    # collapse related contigs and write classification to csv
    collapsed = collapse_contig_classifications(target, classifications, 
                                                csv_writer)

    # convert to recarray and add to master dict
    datatypes = [('contig', '|S32'), ('family', '|S64'), ('type', '|S2')]
    results[target] = np.array(collapsed, dtype=datatypes)

def classify_contigs(contigs, protein_families):
    """Classify contigs based on their domains"""
    # create a list to keep track of classifications
    classifications = []
    
    # classify contigs
    for contig in contigs:
        for protein_family in protein_families:
            if contig.in_family(protein_family):
                classification = [contig.name, protein_family.name, 
                                  protein_family.type]
                # add classification decision to list
                classifications.append(classification)
                
    return classifications

def load_contigs(target, recarray, domains):
    """Loads matching contigs"""
    contigs = []
    
    # keep a record of all the protein domains for each species
    domains[target] = []
    
    for contig_id in set(recarray['target_name']):
        contig_domains = []
        
        # add unique domains associated with the contig id
        for row in recarray[recarray['target_name'] == contig_id]:
            if row['query_name'] not in [d.name for d in contig_domains]:
                domain = ProteinDomain(row['query_name'], row['Evalue'])
                contig_domains.append(domain)
                domains[target].append(row['query_name'])
        
        # create contig and add to the list
        contigs.append(Contig(contig_id, contig_domains))
        
    return contigs

def create_output_dirs():
    """Creates output directories if they don't already exist"""
    for d in ["classifications", "domains", "taps", "taps/gratol", "taps/all"]:
        dir = os.path.join("../csv/", d)
        if not os.path.isdir(dir):
            os.makedirs(dir)

def collapse_contig_classifications(target_name, classifications, csv_writer):
    """Collapses _cN_seqN Trinitiy contig classifications"""
    
    collapsed = []
    
    for classification in classifications:
        # before
        # contig_<Number>_cN_seqN_<translation frame>
        parts = classification[0].split("_")
        general_name = parts[0]
        
        family = classification[1]
        
        # after
        # contig_<Number>_<translation frame>
        new_name = parts[0] + "_" + parts[-1]
        
        # find all similar contigs
        similar = filter(lambda x: x[0].startswith(general_name), 
                         collapsed)
        
        # collapse list of domains for all matches
        matches = set([x[1] for x in similar])
        
        # as long as there are no conflicts, add to new set
        if len(matches) == 0:
            row = tuple([new_name] + classification[1:])
            collapsed.append(row)
            csv_writer.writerow(row)
        elif family in matches:
            # if the same classification already exists, simply ignore for
            # now
            pass
        else:
            print("(%s) multiple classifications for %s: %s" % 
                  (target_name, 
                   general_name, ", ".join(matches.union([family]))))
        
    return collapsed

def get_tap_diversity(results):
    """Generates a dict of TAP diversity organized by target species"""
    # create a dict to keep track unique matches:
    #
    # {
    #    species1: {TR: {...}, TF:{...}, PT: {...}},
    #    species2:...
    #    etc,..
    # }
    #
    tap_diversity = {}
    
    # write summary file with unique and total classification numbers
    for name, cls in results.items():
        # keep track of unique TAP families
        tap_diversity[name] = {
            "TR": set(cls[cls['type'] == "TR"]['family']),
            "TF": set(cls[cls['type'] == "TF"]['family']),
            "PT": set(cls[cls['type'] == "PT"]['family'])
        }
        
    return tap_diversity

def write_classification_summary(results, tap_diversity):
    """Write a summary csv reports"""
    filepath = '../csv/classifications/classification_SUMMARY.csv'
    writer = csv.writer(open(filepath, 'wt'))
    writer.writerow(['species', 'TR (total)', 'TF (total)', 'PT (total)',
                     'TR (unique)', 'TF (unique)', 'PT (unique)'])
    
    # write summary file with unique and total classification numbers
    for name, cls in results.items():
        # add row to summary csv
        writer.writerow((name,
             list(cls['type']).count("TR"), 
             list(cls['type']).count("TF"),
             list(cls['type']).count("PT"),
             len(tap_diversity[name]["TR"]),
             len(tap_diversity[name]["TF"]),
             len(tap_diversity[name]["PT"])
        ))

def write_tap_correlation_csv(taps, prefix='../csv/taps/gratol/correlation', 
                              normalize=False):
    """Generates a CSV file corresponding to the pair-wise correlations between
    the TAP complements of each target species.
    
    Two different calculations are possible depending on whether "normalization"
    is set to True or not:
    
    (1) Normalized (symmetric):
        (A Intersection B) / (A Union B)
    (2) Non-normalized:
        (A Intersection B) / A
        (A Intersection B) / B
    """
    # write additional csv files with pair-wise correlations
    species = taps.keys()

    for category in ["TR", "TF", "PT", "TOTAL"]:
        
        if normalize is True:
            filepath = prefix + ("_normalized_%s.csv" % category)
        else:
            filepath = prefix + ("_%s.csv" % category)

        writer = csv.writer(open(filepath, 'wt'))
        writer.writerow([None] + species)     
        
        for i, species1 in enumerate(species):
            if category == "TOTAL":
                set1 = (taps[species1]['TR'].union(taps[species1]['TF'])
                                            .union(taps[species1]['PT']))
            else:
                set1 = taps[species1][category]
            
            row = [species1]
            
            # find correlation for TAP complements for each species pair
            for j, species2 in enumerate(species):
                # only include one side of diagonal for normalized correlations
                if normalize and j > i:
                    row.append(None)
                    continue

                if category == "TOTAL":
                    set2 = (taps[species2]['TR'].union(taps[species2]['TF'])
                                                .union(taps[species2]['PT']))
                else:
                    set2 = taps[species2][category]
                    
                # calculate union and intersections
                intersection = float(len(set1.intersection(set2)))
                union = float(len(set1.union(set2)))
                
                # if at least some TAPs are shared between species, calculate
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

            writer.writerow(row)

def write_extended_tap_correlation_csv(tap_diversity):
    """Similar to write_tap_correlation_csv, but includes additional TAP
    information from supplement 4 (original paper classifications)."""
    # read in original classifications
    num_species = 20
    
    datatypes = ["S32", "S16"] + (num_species * ['i4'])
    data = np.genfromtxt('../input/Suppl.Table4b.csv', 
                         delimiter=';', names=True, dtype=datatypes)
    
    # get a list of species
    targets = data.dtype.names[2:]
    
    # list of TAP families
    taps = data['TAP']
    
    # add classifications to existing tap_diversity list
    for species in targets:
        tap_diversity[species] = {
            "TR": [],
            "TF": [],
            "PT": []
        }
        
        for i, tap in enumerate(taps):
            # skip no_family_found
            if tap == "no_family_found":
                continue
            
            # otherwise add TAP to set if at least one protein matched
            type_ = data[i][1]
            if data[species][i] > 0:
                tap_diversity[species][type_].append(tap)
                
        # convert to sets
        tap_diversity[species] = {
            "TR": set(tap_diversity[species]['TR']),
            "TF": set(tap_diversity[species]['TF']),
            "PT": set(tap_diversity[species]['PT'])
        }
    
    # write out
    write_tap_correlation_csv(tap_diversity, '../csv/taps/all/correlation_ext')
    write_tap_correlation_csv(tap_diversity, '../csv/taps/all/correlation_ext', 
                              normalize=True)

def write_domain_summary_csv(domains):
    """Writes a sumary of the protein domains matched for each target species"""
    filepath = '../csv/domains/tap_domain_summary.csv'
    writer = csv.writer(open(filepath, 'wt'))
    writer.writerow(["Species", "Domains (total)", "Domains (unique)"])
    
    for species, matches in domains.items():
        writer.writerow([species, len(matches), len(set(matches))])
        
    # pair-wise domain comparison
    species = domains.keys()
    
    filepath = '../csv/domains/tap_domain_correlations.csv'
    writer = csv.writer(open(filepath, 'wt'))
    writer.writerow([None] + species)     
    
    for i, species1 in enumerate(species):
        set1 = set(domains[species1])
        
        row = [species1]
        
        # find correlation for TAP complements for each species pair
        for j, species2 in enumerate(species):
            # only include one side of diagonal
            if j > i:
                row.append(None)
                continue

            set2 = set(domains[species2])
                
            # correlation = INTERSECTION(TAPs) / UNION(TAPs)
            correlation = (len(set1.intersection(set2)) / 
                           float(len(set1.union(set2))))
            row.append(correlation)

        writer.writerow(row)

def write_phylip_char_matrix(targets, protein_families, outgroup=None, 
                             tap_classes=None, filename="infile"):
    """Writes a phylip boolean character matrix representing the presence or
    absense of each TAP family for each species."""
    # include all TAP classes by default
    if tap_classes is None:
        tap_classes = ["TF", "TR", "PT"]
        
    # make a copy of species dictionary so that when outgroup is removed it
    # does not affect original dict
    species = targets.copy()
    
    families = [x.name for x in protein_families if x.type in tap_classes]
    
    # open file for output
    fp = open(os.path.join("../trees/", filename), "w")
    fp.write("%d %d\n" % (len(species), len(families)))
    
    # put outgroup first
    if outgroup is not None:
        # add outgroup
        row = get_char_matrix_row(outgroup, species[outgroup], families, tap_classes)
        fp.write(outgroup.ljust(28) + "".join(row) + "\n")
        
        # remove from set so it doesn't get double-counted
        species.pop(outgroup)
    
    # process rest of the species
    for name, taps in species.items():
        row = get_char_matrix_row(name, taps, families, tap_classes)
                
        fp.write(name.ljust(28) + "".join(row) + "\n")
        
def get_char_matrix_row(species, taps, families, tap_classes):
    """Computes a single row in the character matrix"""
    row = []
    
    matched_families = set([])
    
    for cls in tap_classes:
        matched_families = matched_families.union(taps[cls])        
    
    # Output a boolean character row
    for family in families:
        if family in matched_families:
            row.append("1")
        else:
            row.append("0")
            
    return row
           
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
        ProteinFamily("BBR/BPC", "TF", ["GAGA_bind"]),
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
        ProteinFamily("DUF246 domain containing", "PT", ["O-FucT"]),
        ProteinFamily("DUF296 domain containing", "PT", ["DUF296"]),
        ProteinFamily("DUF547 domain containing", "PT", ["DUF547"]),
        ProteinFamily("DUF632 domain containing", "PT", ["DUF632"]),
        ProteinFamily("DUF833 domain containing", "PT", ["NRDE"]),
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
        ProteinFamily("IWS1", "TR", ["Med26"], []),
        ProteinFamily("Jumonji", "TR", ["JmjC", "JmjN"], ["ARID", "GATA", "zf-C2H2", "Alfin-like"]),
        ProteinFamily("LFY", "TF", ["FLO_LFY"]),
        ProteinFamily("LIM", "TF", ["two_or_more_LIM"]),
        ProteinFamily("LUG", "TR", ["LUFS_Domain"]),
        ProteinFamily("MADS", "TF", ["SRF-TF"]),
        ProteinFamily("MBF1", "TR", ["MBF1"]),
        ProteinFamily("MED6", "TR", ["Med6"]),
        ProteinFamily("MED7", "TR", ["Med7"]),
        ProteinFamily("mTERF", "TF", ["mTERF"]),
        ProteinFamily("MYB", "TF", ["two_or_more_Myb_DNA-binding"], ["ARID", "G2-like_Domain", "Response_reg", "trihelix"]),
        ProteinFamily("MYB-related", "TF", ["Myb_DNA-binding"], ["ARID", "G2-like_Domain", "Response_reg", "trihelix", "two_or_more_Myb_DNA-binding"]),
        ProteinFamily("NAC", "TF", ["NAC_plant"]),
        ProteinFamily("NZZ", "TF", ["NOZZLE"]),
        ProteinFamily("OFP", "TR", ["Ovate"]),
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
        ProteinFamily("SOH1", "TR", ["Med31"]),
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

    def in_family(self, family):
        """Checks to see whether the contig belongs to a protein family.
        
        Parameters:
        -----------
        family : ProteinFamily
            Protein family classification rules
        """
        contig_domains = self.get_domain_names()
        
        return (family.requires.issubset(contig_domains) and 
                family.forbids.isdisjoint(contig_domains) and
                (len(family.alt_domains) == 0 or 
                 len(family.alt_domains.intersection(contig_domains)) > 0))
    
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
    
if __name__ == '__main__':
    #sys.exit(main())
    results = main()


