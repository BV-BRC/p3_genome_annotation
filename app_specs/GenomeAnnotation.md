
# Application specification: GenomeAnnotation

This is the application specification for service with identifier GenomeAnnotation.

The backend script implementing the application is [App-GenomeAnnotation](../service-scripts/App-GenomeAnnotation.pm).

The raw JSON file for this specification is [GenomeAnnotation.json](GenomeAnnotation.json).

This service performs the following task:   Calls genes and functionally annotate input contig set.

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | :------: | ------------ |
| contigs | Contig file | WS: Contigs  | :heavy_check_mark: |  |
| scientific_name | Scientific Name | string  | :heavy_check_mark: |  |
| taxonomy_id | NCBI Taxonomy ID | int  |  |  |
| code | Genetic Code | int  | :heavy_check_mark: | 0 |
| domain | Domain | enum  | :heavy_check_mark: | auto |
| public | Public | bool  |  | 0 |
| queue_nowait | Don't wait on indexing queue | bool  |  | 0 |
| skip_indexing | Don't index genome | bool  |  | 0 |
| skip_workspace_output | Don't write to workspace | bool  |  | 0 |
| output_path | Output Folder | folder  |  |  |
| output_file | File Basename | wsid  |  |  |
| reference_genome_id | Reference genome ID | string  |  |  |
| reference_virus_name | Reference virus name | string  |  |  |
| container_id | (Internal) Container to use for this run | string  |  |  |
| indexing_url | (Internal) Override Data API URL for use in indexing | string  |  |  |
| _parent_job | (Internal) Parent job for this annotation | string  |  |  |
| fix_errors | Automatically fix errors? | bool  |  |  |
| fix_frameshifts | Fix frameshifts? | bool  |  |  |
| enable_debug | Enable debug? | bool  |  |  |
| verbose_level | Set verbose level | int  |  |  |
| workflow | Custom workflow | string  |  |  |
| recipe | Annotation recipe | string  |  |  |
| disable_replication | Disable replication? | bool  |  |  |
| analyze_quality | Enable quality analysis of genome | bool  |  |  |
| assembly_output | Workspace path holding assembly output for this genome | folder  |  |  |
| custom_pipeline | Customize the RASTtk pipeline | group  |  |  |

