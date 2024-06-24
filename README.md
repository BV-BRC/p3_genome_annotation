# Genome Annotation

This module defines the [BV-BRC](https://bv-brc.org) service script for the genome annotation services.

There are three distinct application types defined for it:

1. [GenomeAnnotation](app_scripts/GenomeAnnotation.md): Service that provides the backend for the BV-BRC web interface; it takes contigs as input.
2. [GenomeAnnotationGenbank](app_scripts/GenomeAnnotationGenbank.md): Currently only available via the `p3-submit-genome-annotation` command-line script; allows the processing of Genbank files as input.
3. [ComprehensiveGenomeAnalysis](app_scripts/ComprehensiveGenomeAnalysis.md): Provides the backend for the BV-BRC Comprehensive Genome Annotation service.


