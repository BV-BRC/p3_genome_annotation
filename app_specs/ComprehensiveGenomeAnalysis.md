
# Application specfication: ComprehensiveGenomeAnalysis

This is the application specification for service with identifier ComprehensiveGenomeAnalysis.

The backend script implementing the application is [App-ComprehensiveGenomeAnalysis](service-scripts/App-ComprehensiveGenomeAnalysis.pm).

This service performs the following task:   Analyze a genome from reads or contigs, generating a detailed analysis report.

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | :------: | ------------ |
|  |  |   | :x: |  |
| input_type | Input Type | enum  | :heavy_check_mark: |  |
| output_path | Output Folder | folder  | :heavy_check_mark: |  |
| output_file | File Basename | wsid  | :heavy_check_mark: |  |
|  |  |   | :x: |  |
| paired_end_libs |  | group  | :x: |  |
| single_end_libs |  | group  | :x: |  |
| srr_ids | SRR ID | string  | :x: |  |
| reference_assembly | Contig file | WS: Contigs  | :x: |  |
| recipe | Assembly recipe | enum  | :x: | auto |
| racon_iter | Racon iterations | int  | :x: | 2 |
| pilon_iter | Pilon iterations | int  | :x: | 2 |
| trim | trim_reads | boolean  | :x: | 0 |
| min_contig_len | Minimal output contig length | int  | :x: | 300 |
| min_contig_cov | Minimal output contig coverage | float  | :x: | 5 |
| genome_size | Genome Size | string  | :x: | 5M |
|  |  |   | :x: |  |
| gto | Preannotated genome object | WS: Genome  | :x: |  |
| genbank_file | Genbank file | WS: genbank_file  | :x: |  |
| contigs | Contig file | WS: Contigs  | :x: |  |
| scientific_name | Scientific Name | string  | :heavy_check_mark: |  |
| taxonomy_id | NCBI Taxonomy ID | int  | :x: |  |
| code | Genetic Code | int  | :heavy_check_mark: | 0 |
| domain | Domain | enum  | :heavy_check_mark: | auto |
|  |  |   | :x: |  |
| public | Public | bool  | :x: | 0 |
| queue_nowait | Don't wait on indexing queue | bool  | :x: | 0 |
| skip_indexing | Don't index genome | bool  | :x: | 0 |
| reference_genome_id | Reference genome ID | string  | :x: |  |
| _parent_job | (Internal) Parent job for this annotation | string  | :x: |  |
| workflow | Custom workflow | string  | :x: |  |
| analyze_quality | Enable quality analysis of genome | bool  | :x: |  |
| debug_level | Debug level | int  | :x: | 0 |

