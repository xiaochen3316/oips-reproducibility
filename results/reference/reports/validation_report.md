# Validation report

Overall status: **PASS**

Pass: 23  
Fail: 0  
Warning: 0

## Checks

| Check ID | Status | Summary | Evidence |
|---|---|---|---|
| checksums | PASS | validated | {"files":31} |
| release_tables | PASS | validated | {"tables":9} |
| data_manifest | PASS | validated | {"rows":30} |
| scientific_coverage | PASS | validated | {"records":1742,"targets":21} |
| expected_summary | PASS | validated | {"keys":24} |
| sensitive_content | PASS | validated | {"classes":[],"files":[],"match_count":0} |
| bundle_contract | PASS | stage contracts validated | {} |
| snapshot_cluster_summary | PASS | cluster snapshot matches | {"mismatches":[]} |
| snapshot_cluster_numeric | PASS | snapshot values match | {"metrics":2} |
| snapshot_cluster_distributions | PASS | snapshot values match | {} |
| snapshot_exemplar_static | PASS | snapshot values match | {"metrics":9} |
| snapshot_exemplar_state | PASS | snapshot values match | {} |
| snapshot_analysis_metrics | PASS | snapshot values match | {"metrics":11} |
| snapshot_bootstrap_family | PASS | snapshot values match | {"metrics":14} |
| snapshot_orel_ablation | PASS | snapshot values match | {"metrics":5} |
| snapshot_supplementary_analyses | PASS | snapshot values match | {"metrics":46} |
| snapshot_single_tool_complete_cases | PASS | snapshot values match | {} |
| snapshot_posthoc_counts | PASS | snapshot values match | {} |
| snapshot_posthoc_distributions | PASS | snapshot values match | {} |
| snapshot_posthoc_reference_case | PASS | snapshot values match | {} |
| snapshot_analysis_counts | PASS | snapshot values match | {} |
| snapshot_representative_order | PASS | representative case order matches | {"actual":["5J89","5TBM","4W9H"]} |
| report_sensitive_content | PASS | no sensitive content detected | {"classes":[],"files":[],"match_count":0} |
