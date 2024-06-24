
# Application specfication: GenomeAnnotation

This is the application specification for service with identifier GenomeAnnotation.

The backend script implementing the application is [App-GenomeAnnotation](service-scripts/App-GenomeAnnotation.pm).

This service performs the following task:   Calls genes and functionally annotate input contig set.

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | :------: | ------------ |
| contigs | Contig file | WS: Contigs  | :heavy_check_mark: |  |
| scientific_name | Scientific Name | string  | :heavy_check_mark: |  |
| taxonomy_id | NCBI Taxonomy ID | int  | :x: |  |
| code | Genetic Code | int  | :heavy_check_mark: | 0 |
| domain | Domain | enum  | :heavy_check_mark: | auto |
| public | Public | bool  | :x: | 0 |
| queue_nowait | Don't wait on indexing queue | bool  | :x: | 0 |
| skip_indexing | Don't index genome | bool  | :x: | 0 |
| skip_workspace_output | Don't write to workspace | bool  | :x: | 0 |
| output_path | Output Folder | folder  | :x: |  |
| output_file | File Basename | wsid  | :x: |  |
| reference_genome_id | Reference genome ID | string  | :x: |  |
| reference_virus_name | Reference virus name | string  | :x: |  |
| container_id | (Internal) Container to use for this run | string  | :x: |  |
| indexing_url | (Internal) Override Data API URL for use in indexing | string  | :x: |  |
| _parent_job | (Internal) Parent job for this annotation | string  | :x: |  |
| fix_errors | Automatically fix errors? | bool  | :x: |  |
| fix_frameshifts | Fix frameshifts? | bool  | :x: |  |
| enable_debug | Enable debug? | bool  | :x: |  |
| verbose_level | Set verbose level | int  | :x: |  |
| workflow | Custom workflow | string  | :x: |  |
| recipe | Annotation recipe | string  | :x: |  |
| disable_replication | Disable replication? | bool  | :x: |  |
| analyze_quality | Enable quality analysis of genome | bool  | :x: |  |
| assembly_output | Workspace path holding assembly output for this genome | folder  | :x: |  |
| custom_pipeline | Customize the RASTtk pipeline | group  | :x: |  |

