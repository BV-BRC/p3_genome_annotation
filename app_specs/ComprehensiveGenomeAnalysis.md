
# Application specfication: ComprehensiveGenomeAnalysis

This is the application specification for service with identifier ComprehensiveGenomeAnalysis.

The backend script implementing the application is [App-ComprehensiveGenomeAnalysis](service-scripts/App-ComprehensiveGenomeAnalysis.pm).

This service performs the following task:   Analyze a genome from reads or contigs, generating a detailed analysis report.

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | :------: | ------------ |
|  |  |   |  |  |
| input_type | Input Type | enum  | :heavy_check_mark: |  |
| output_path | Output Folder | folder  | :heavy_check_mark: |  |
| output_file | File Basename | wsid  | :heavy_check_mark: |  |
|  |  |   |  |  |
| paired_end_libs |  | group  |  |  |
| single_end_libs |  | group  |  |  |
| srr_ids | SRR ID | string  |  |  |
| reference_assembly | Contig file | WS: Contigs  |  |  |
| recipe | Assembly recipe | enum  |  | auto |
| racon_iter | Racon iterations | int  |  | 2 |
| pilon_iter | Pilon iterations | int  |  | 2 |
| trim | trim_reads | boolean  |  | 0 |
| min_contig_len | Minimal output contig length | int  |  | 300 |
| min_contig_cov | Minimal output contig coverage | float  |  | 5 |
| genome_size | Genome Size | string  |  | 5M |
|  |  |   |  |  |
| gto | Preannotated genome object | WS: Genome  |  |  |
| genbank_file | Genbank file | WS: genbank_file  |  |  |
| contigs | Contig file | WS: Contigs  |  |  |
| scientific_name | Scientific Name | string  | :heavy_check_mark: |  |
| taxonomy_id | NCBI Taxonomy ID | int  |  |  |
| code | Genetic Code | int  | :heavy_check_mark: | 0 |
| domain | Domain | enum  | :heavy_check_mark: | auto |
|  |  |   |  |  |
| public | Public | bool  |  | 0 |
| queue_nowait | Don't wait on indexing queue | bool  |  | 0 |
| skip_indexing | Don't index genome | bool  |  | 0 |
| reference_genome_id | Reference genome ID | string  |  |  |
| _parent_job | (Internal) Parent job for this annotation | string  |  |  |
| workflow | Custom workflow | string  |  |  |
| analyze_quality | Enable quality analysis of genome | bool  |  |  |
| debug_level | Debug level | int  |  | 0 |

