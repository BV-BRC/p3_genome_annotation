
# Application specfication: GenomeAnnotationGenbankTest

This is the application specification for service with identifier GenomeAnnotationGenbankTest.

The backend script implementing the application is [App-GenomeAnnotationGenbankTest](service-scripts/App-GenomeAnnotationGenbankTest.pm).

This service performs the following task:   Calls genes and functionally annotate input contig set.

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | :------: | ------------ |
| genbank_file | Genbank file | WS: genbank_file  | :heavy_check_mark: |  |
| public | Public | bool  | :x: | 0 |
| queue_nowait | Don't wait on indexing queue | bool  | :x: | 0 |
| output_path | Output Folder | folder  | :x: |  |
| output_file | File Basename | wsid  | :x: |  |
| fix_errors | Automatically fix errors? | bool  | :x: |  |
| fix_frameshifts | Fix frameshifts? | bool  | :x: |  |
| enable_debug | Enable debug? | bool  | :x: |  |
| verbose_level | Set verbose level | int  | :x: |  |
| disable_replication | Disable replication? | bool  | :x: |  |
| custom_pipeline | Customize the RASTtk pipeline | group  | :x: |  |

