# Genome Annotation

## Overview

Genome annotation is the process of identifying functional elements
along the sequence of a genome. The Genome Annotation Service in
BV-BRC uses the RAST tool kit (RASTtk) to annotate genomic features
in bacteria, the Viral Genome ORF Reader (VIGOR4) to annotate
viruses, and PHANOTATE to annotate bacteriophages. All public
genomes in BV-BRC have been consistently annotated with these
services. Researchers can also submit their own private genome to the
annotation service, where it will be deposited into their private
workspace for their perusal and comparison with other BV-BRC genomes
using BV-BRC analysis tools and services, e.g., Phylogenetic Tree,
Genome Alignment, Protein Family Sorter, Comparative Pathway Viewer,
Similar Genome Finder.

This module defines the [BV-BRC](https://bv-brc.org) service script for the genome annotation services.

There are three distinct application types defined here:

1. [GenomeAnnotation](app_specs/GenomeAnnotation.md): Service that provides the backend for the BV-BRC web interface; it takes contigs as input.
2. [GenomeAnnotationGenbank](app_specs/GenomeAnnotationGenbank.md): Currently only available via the `p3-submit-genome-annotation` command-line script; allows the processing of Genbank files as input.
3. [ComprehensiveGenomeAnalysis](app_specs/ComprehensiveGenomeAnalysis.md): Provides the backend for the BV-BRC Comprehensive Genome Annotation service.

## About this module

This module is a component of the BV-BRC build system. It is designed to fit into the
`dev_container` infrastructure which manages development and production deployment of
the components of the BV-BRC. More documentation is available [here](https://github.com/BV-BRC/dev_container/tree/master/README.md).

The code in this module provides the BV-BRC application service wrapper scripts for the genome annotation service as well
as some backend utilities:

| Script name | Purpose |
| ----------- | ------- |
| [App-ComprehensiveGenomeAnnotation.pl](service-scripts/App-ComprehensiveGenomeAnnotation.pl) | App script for the comprehnsive genome annotation service |

## See also

* [Genome Annotation Service](https://bv-brc.org/app/Annotation)
* [Quick Reference](https://www.bv-brc.org/docs/quick_references/services/genome_annotation_service.html)
* [Genome Annotation Service Tutorial](https://www.bv-brc.org/docs/tutorial/genome_annotation/genome_annotation.html)



## References

Aziz, R. K. et al. The RAST Server: rapid annotations using subsystems technology. BMC genomics 9, 75 (2008).

VIGOR4, https://github.com/JCVenterInstitute/VIGOR4.

McNair, K. et al. in Bacteriophages 231-238 (Springer, 2018).

McNair, K., Zhou, C., Dinsdale, E. A., Souza, B. & Edwards, R. A. PHANOTATE: a novel approach to gene identification in phage genomes. Bioinformatics 35, 4537-4542 (2019).

Overbeek, R. et al. The subsystems approach to genome annotation and its use in the project to annotate 1000 genomes. Nucleic acids research 33, 5691-5702 (2005).

Overbeek, R. et al. The SEED and the Rapid Annotation of microbial genomes using Subsystems Technology (RAST). 42, D206-D214 (2013).

Davis, J. J. et al. The PATRIC Bioinformatics Resource Center: expanding data and analysis capabilities. Nucleic acids research 48, D606-D612 (2020).

Brettin, T. et al. RASTtk: a modular and extensible implementation of the RAST algorithm for building custom annotation pipelines and annotating batches of genomes. Scientific reports 5, 8365 (2015).

Mao, C. et al. Curation, integration and visualization of bacterial virulence factors in PATRIC. Bioinformatics 31, 252-258 (2015).

Ye, J., McGinnis, S. & Madden, T. L. BLAST: improvements for better sequence analysis. Nucleic acids research 34, W6-W9 (2006).

Croucher, N. J., Vernikos, G. S., Parkhill, J. & Bentley, S. D. Identification, variation and transcription of pneumococcal repeat sequences. BMC genomics 12, 1-13 (2011).

Hyatt, D. et al. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC bioinformatics 11, 1-11 (2010).

Delcher, A. L., Bratke, K. A., Powers, E. C. & Salzberg, S. L. Identifying bacterial genes and endosymbiont DNA with Glimmer. Bioinformatics 23, 673-679 (2007).

Davis, J. J. et al. Antimicrobial resistance prediction in PATRIC and RAST. Scientific reports 6, 27930 (2016).

Kent, W. J. BLAT—the BLAST-like alignment tool. Genome research 12, 656-664 (2002).

Johnson, M. et al. NCBI BLAST: a better web interface. Nucleic acids research 36, W5-W9 (2008).

Liu, B., Zheng, D., Jin, Q., Chen, L. & Yang, J. VFDB 2019: a comparative pathogenomic platform with an interactive web interface. Nucleic acids research 47, D687-D692 (2019).

Xiang, Z. et al. VIOLIN: vaccine investigation and online information network. Nucleic acids research 36, D923-D928 (2007).

Alcock, B. P. et al. CARD 2020: antibiotic resistome surveillance with the comprehensive antibiotic resistance database. Nucleic acids research 48, D517-D525 (2020).

Liu, B. & Pop, M. ARDB—antibiotic resistance genes database. Nucleic acids research 37, D443-D447 (2009).

Antonopoulos, D. A. et al. PATRIC as a unique resource for studying antimicrobial resistance. Briefings in bioinformatics (2017).

Saier Jr, M. H. et al. The transporter classification database (TCDB): recent advances. Nucleic acids research 44, D372-D379 (2016).

Wishart, D. S. et al. DrugBank 5.0: a major update to the DrugBank database for 2018. Nucleic acids research 46, D1074-D1082 (2018).

Chen, X., Ji, Z. L. & Chen, Y. Z. TTD: therapeutic target database. Nucleic acids research 30, 412-415 (2002).

Davis, J. J. et al. PATtyFams: Protein families for the microbial genomes in the PATRIC database. 7, 118 (2016).

Akhter, S., Aziz, R. K. & Edwards, R. A. PhiSpy: a novel algorithm for finding prophages in bacterial genomes that combines similarity-and composition-based strategies. Nucleic acids research 40, e126-e126 (2012).

Osawa, S., Jukes, T. H., Watanabe, K. & Muto, A. Recent evidence for evolution of the genetic code. Microbiological reviews 56, 229-264 (1992).

Rivest, R. & Dusse, S. (MIT Laboratory for Computer Science Cambridge, 1992).

Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P. & Tyson, G. W. CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome research 25, 1043-1055 (2015).