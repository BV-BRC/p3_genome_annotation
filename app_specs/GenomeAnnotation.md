
# Application specficiation: GenomeAnnotation

This is the application specification for service with identifier GenomeAnnotation.

The backend script implementing the application is [App-GenomeAnnotation](service-scripts/App-GenomeAnnotation.pm).

This service does the following:

   Calls genes and functionally annotate input contig set.

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | -------- | ------------ |
| contigs | Contig file | WS: Contigs  | YES |  |
| scientific_name | Scientific Name | string  | YES |  |
| taxonomy_id | NCBI Taxonomy ID | int  | NO |  |
| code | Genetic Code | int  | YES | 0 |
| domain | Domain | enum  | YES | auto |
| public | Public | bool  | NO | 0 |
| queue_nowait | Don't wait on indexing queue | bool  | NO | 0 |
| skip_indexing | Don't index genome | bool  | NO | 0 |
| skip_workspace_output | Don't write to workspace | bool  | NO | 0 |
| output_path | Output Folder | folder  | NO |  |
| output_file | File Basename | wsid  | NO |  |
| reference_genome_id | Reference genome ID | string  | NO |  |
| reference_virus_name | Reference virus name | string  | NO |  |
| container_id | (Internal) Container to use for this run | string  | NO |  |
| indexing_url | (Internal) Override Data API URL for use in indexing | string  | NO |  |
| _parent_job | (Internal) Parent job for this annotation | string  | NO |  |
| fix_errors | Automatically fix errors? | bool  | NO |  |
| fix_frameshifts | Fix frameshifts? | bool  | NO |  |
| enable_debug | Enable debug? | bool  | NO |  |
| verbose_level | Set verbose level | int  | NO |  |
| workflow | Custom workflow | string  | NO |  |
| recipe | Annotation recipe | string  | NO |  |
| disable_replication | Disable replication? | bool  | NO |  |
| analyze_quality | Enable quality analysis of genome | bool  | NO |  |
| assembly_output | Workspace path holding assembly output for this genome | folder  | NO |  |
| custom_pipeline | Customize the RASTtk pipeline | group  | NO |  |

