Rotation Project: Green Algae TAP Family Evolution
==================================================

Overview
--------
This repository includes all of the code related to a project I worked on as 
part of a first-year lab rotation in the [Delwiche laboratory](http://www.life.umd.edu/labs/delwiche/home.html).

![Spirogyra](https://raw.github.com/khughitt/rotation_delwiche/master/extra/3x2_millimeters_of_Spirogyra.jpg)
(source: [Wikipedia](http://commons.wikimedia.org/wiki/File:3x2_millimeters_of_Spirogyra.jpg))

The project was inspired by a [2009 study by Lang et. al.](http://www.ncbi.nlm.nih.gov/pubmed/20644220)
that, among other things,  looked at correlations between increases in the 
Transcription-Associated Protein (TAP) family diversity and several milestones 
in plant evolution such as the development of land plants and the radiation 
of angiosperms. While the Lang et al. study included genome data from a large 
number  of plant species, there were several lineages for which no data was 
available, including the [Charophyte green algae](http://en.wikipedia.org/wiki/Charophyta).

Species studied:

* Cephaleuors parasiticus 
* Mesostigma viride
* Nitella mirabilis
* Phaeophila dendroides
* Pyramimonas parkeae
* Spirogyra pratensis

The goal for this project was to extend the analyses done by Lang et al.
to include additional [green algae](http://en.wikipedia.org/wiki/Green_algae) 
datasets, including several Charophyta species.

Use of RNA-Seq Data
-------------------
An additional aim of the project was to determine the feasibility of using
expression ([RNA-Seq](http://en.wikipedia.org/wiki/RNA-Seq)) data to measure
TAP diversity.

The original study included only species for which complete genome sequences
were available. For this project, we had not genome sequences, but large 
transcriptome data-sets. The obviously has some major implications for comparing
the results of this project with results from the previous study, and also
comparing different species within the project.

The result is that we can only make claims about the *minimum* number of TAP 
families present in a given species.

RNA-Seq data used for this project was assembled de novo using [Trinity](http://trinityrnaseq.sourceforge.net/).

Finding TAP domains
-------------------
The first step in the project was to find RNA-Seq contigs containing protein 
domains found in TAPs.

[Profile Hidden Markov Models (HMMs)](http://www.ncbi.nlm.nih.gov/pubmed/9918945)
were  downloaded for each of the TAP domains listed in the Lang et al. study.
Many of the profile HMMs were downloaded from the [Pfam](http://pfam.sanger.ac.uk/)
database, along with several custom profile HMMs created by Lang et al.

In the several years that have elapsed since the original study, a number of
Pfam domains identifiers have changed, and it was necessary to map the ids
used in the original study to the new identifiers:

    DUF1004 -> GAGA_bind
    DUF246 -> O-FucT
    DUF623 -> Ovate
    DUF833 -> NRDE
    IWS1_C -> Med26
    MED6 -> Med6
    MED7 -> Med7
    SOH1 -> Med31
    
Old names are listed on the left side while the newer Pfam identifiers are on
the right side.

HMMER](http://hmmer.janelia.org/) was used to scan for protein domains in the
RNA-Seq contig data.

Todo
----
* Decide on how to deal with multiple-classification contigs
