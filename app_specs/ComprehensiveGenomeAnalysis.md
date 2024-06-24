
# Application specficiation: ComprehensiveGenomeAnalysis

This is the application specification for service with identifier ComprehensiveGenomeAnalysis.

The backend script implementing the application is [App-ComprehensiveGenomeAnalysis](service-scripts/App-ComprehensiveGenomeAnalysis.pm).

This service does the following:

   Analyze a genome from reads or contigs, generating a detailed analysis report.

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | -------- | ------------ |
|  |  |   | NO |  |
| input_type | Input Type | enum  | YES |  |
| output_path | Output Folder | folder  | YES |  |
| output_file | File Basename | wsid  | YES |  |
|  |  |   | NO |  |
| paired_end_libs |  | group  | NO |  |
| single_end_libs |  | group  | NO |  |
| srr_ids | SRR ID | string  | NO |  |
| reference_assembly | Contig file | WS: Contigs  | NO |  |
| recipe | Assembly recipe | enum  | NO | auto |
| racon_iter | Racon iterations | int  | NO | 2 |
| pilon_iter | Pilon iterations | int  | NO | 2 |
| trim | trim_reads | boolean  | NO | 0 |
| min_contig_len | Minimal output contig length | int  | NO | 300 |
| min_contig_cov | Minimal output contig coverage | float  | NO | 5 |
| genome_size | Genome Size | string  | NO | 5M |
|  |  |   | NO |  |
| gto | Preannotated genome object | WS: Genome  | NO |  |
| genbank_file | Genbank file | WS: genbank_file  | NO |  |
| contigs | Contig file | WS: Contigs  | NO |  |
| scientific_name | Scientific Name | string  | YES |  |
| taxonomy_id | NCBI Taxonomy ID | int  | NO |  |
| code | Genetic Code | int  | YES | 0 |
| domain | Domain | enum  | YES | auto |
|  |  |   | NO |  |
| public | Public | bool  | NO | 0 |
| queue_nowait | Don't wait on indexing queue | bool  | NO | 0 |
| skip_indexing | Don't index genome | bool  | NO | 0 |
| reference_genome_id | Reference genome ID | string  | NO |  |
| _parent_job | (Internal) Parent job for this annotation | string  | NO |  |
| workflow | Custom workflow | string  | NO |  |
| analyze_quality | Enable quality analysis of genome | bool  | NO |  |
| debug_level | Debug level | int  | NO | 0 |

