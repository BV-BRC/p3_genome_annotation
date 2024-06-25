
# Application specification: GenomeAnnotationGenbank

This is the application specification for service with identifier GenomeAnnotationGenbank.

The backend script implementing the application is [App-GenomeAnnotationGenbank.pl](../service-scripts/App-GenomeAnnotationGenbank.pl).

The raw JSON file for this specification is [GenomeAnnotationGenbank.json](GenomeAnnotationGenbank.json).

This service performs the following task:   Calls genes and functionally annotate input contig set.

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | :------: | ------------ |
| genbank_file | Genbank file | WS: genbank_file  | :heavy_check_mark: |  |
| public | Public | bool  |  | 0 |
| queue_nowait | Don't wait on indexing queue | bool  |  | 0 |
| skip_indexing | Don't index genome | bool  |  | 0 |
| skip_workspace_output | Don't write to workspace | bool  |  | 0 |
| container_id | (Internal) Container to use for this run | string  |  |  |
| indexing_url | (Internal) Override Data API URL for use in indexing | string  |  |  |
| output_path | Output Folder | folder  |  |  |
| output_file | File Basename | wsid  |  |  |
| reference_virus_name | Reference virus name | string  |  |  |
| workflow | Custom workflow | string  |  |  |
| recipe | Annotation recipe | string  |  |  |
| scientific_name | Scientific Name | string  |  |  |
| taxonomy_id | NCBI Taxonomy ID | int  |  |  |
| code | Genetic Code | enum  |  |  |
| domain | Domain | enum  |  | Bacteria |
| import_only | Import only | bool  |  |  |
| raw_import_only | Raw import only | bool  |  |  |
| skip_contigs | Skip contigs on import | bool  |  |  |
| fix_errors | Automatically fix errors? | bool  |  |  |
| fix_frameshifts | Fix frameshifts? | bool  |  |  |
| enable_debug | Enable debug? | bool  |  |  |
| verbose_level | Set verbose level | int  |  |  |
| disable_replication | Disable replication? | bool  |  |  |
| custom_pipeline | Customize the RASTtk pipeline | group  |  |  |

