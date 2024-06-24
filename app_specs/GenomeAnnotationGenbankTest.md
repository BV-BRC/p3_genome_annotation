
# Application specfication: GenomeAnnotationGenbankTest

This is the application specification for service with identifier GenomeAnnotationGenbankTest.

The backend script implementing the application is [App-GenomeAnnotationGenbankTest](service-scripts/App-GenomeAnnotationGenbankTest.pm).

This service performs the following task:   Calls genes and functionally annotate input contig set.

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | -------- | ------------ |
| genbank_file | Genbank file | WS: genbank_file  | YES |  |
| public | Public | bool  | NO | 0 |
| queue_nowait | Don't wait on indexing queue | bool  | NO | 0 |
| output_path | Output Folder | folder  | NO |  |
| output_file | File Basename | wsid  | NO |  |
| fix_errors | Automatically fix errors? | bool  | NO |  |
| fix_frameshifts | Fix frameshifts? | bool  | NO |  |
| enable_debug | Enable debug? | bool  | NO |  |
| verbose_level | Set verbose level | int  | NO |  |
| disable_replication | Disable replication? | bool  | NO |  |
| custom_pipeline | Customize the RASTtk pipeline | group  | NO |  |

