
# Application specfication: GenomeAnnotationGenbank

This is the application specification for service with identifier GenomeAnnotationGenbank.

The backend script implementing the application is [App-GenomeAnnotationGenbank](service-scripts/App-GenomeAnnotationGenbank.pm).

This service performs the following task:   Calls genes and functionally annotate input contig set.

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | :------: | ------------ |
| genbank_file | Genbank file | WS: genbank_file  | :heavy_check_mark: |  |
| public | Public | bool  | :x: | 0 |
| queue_nowait | Don't wait on indexing queue | bool  | :x: | 0 |
| skip_indexing | Don't index genome | bool  | :x: | 0 |
| skip_workspace_output | Don't write to workspace | bool  | :x: | 0 |
| container_id | (Internal) Container to use for this run | string  | :x: |  |
| indexing_url | (Internal) Override Data API URL for use in indexing | string  | :x: |  |
| output_path | Output Folder | folder  | :x: |  |
| output_file | File Basename | wsid  | :x: |  |
| reference_virus_name | Reference virus name | string  | :x: |  |
| workflow | Custom workflow | string  | :x: |  |
| recipe | Annotation recipe | string  | :x: |  |
| scientific_name | Scientific Name | string  | :x: |  |
| taxonomy_id | NCBI Taxonomy ID | int  | :x: |  |
| code | Genetic Code | enum  | :x: |  |
| domain | Domain | enum  | :x: | Bacteria |
| import_only | Import only | bool  | :x: |  |
| raw_import_only | Raw import only | bool  | :x: |  |
| skip_contigs | Skip contigs on import | bool  | :x: |  |
| fix_errors | Automatically fix errors? | bool  | :x: |  |
| fix_frameshifts | Fix frameshifts? | bool  | :x: |  |
| enable_debug | Enable debug? | bool  | :x: |  |
| verbose_level | Set verbose level | int  | :x: |  |
| disable_replication | Disable replication? | bool  | :x: |  |
| custom_pipeline | Customize the RASTtk pipeline | group  | :x: |  |

