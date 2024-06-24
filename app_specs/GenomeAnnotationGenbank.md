
# Application specficiation: GenomeAnnotationGenbank

This is the application specification for service with identifier GenomeAnnotationGenbank.

The backend script implementing the application is [App-GenomeAnnotationGenbank](service-scripts/App-GenomeAnnotationGenbank.pm).

This service does the following:

   Calls genes and functionally annotate input contig set.

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | -------- | ------------ |
| genbank_file | Genbank file | WS: genbank_file  | YES |  |
| public | Public | bool  | NO | 0 |
| queue_nowait | Don't wait on indexing queue | bool  | NO | 0 |
| skip_indexing | Don't index genome | bool  | NO | 0 |
| skip_workspace_output | Don't write to workspace | bool  | NO | 0 |
| container_id | (Internal) Container to use for this run | string  | NO |  |
| indexing_url | (Internal) Override Data API URL for use in indexing | string  | NO |  |
| output_path | Output Folder | folder  | NO |  |
| output_file | File Basename | wsid  | NO |  |
| reference_virus_name | Reference virus name | string  | NO |  |
| workflow | Custom workflow | string  | NO |  |
| recipe | Annotation recipe | string  | NO |  |
| scientific_name | Scientific Name | string  | NO |  |
| taxonomy_id | NCBI Taxonomy ID | int  | NO |  |
| code | Genetic Code | enum  | NO |  |
| domain | Domain | enum  | NO | Bacteria |
| import_only | Import only | bool  | NO |  |
| raw_import_only | Raw import only | bool  | NO |  |
| skip_contigs | Skip contigs on import | bool  | NO |  |
| fix_errors | Automatically fix errors? | bool  | NO |  |
| fix_frameshifts | Fix frameshifts? | bool  | NO |  |
| enable_debug | Enable debug? | bool  | NO |  |
| verbose_level | Set verbose level | int  | NO |  |
| disable_replication | Disable replication? | bool  | NO |  |
| custom_pipeline | Customize the RASTtk pipeline | group  | NO |  |

