# Public OIPS data snapshot

This directory contains the reviewed 1,742-record, 21-target public input snapshot for the OIPS reproducibility workflow.

`static/tool_pocket_features.csv` is the allow-listed detector feature table. `metadata/systems.tsv` records approved project annotations and leaves unverified preparation or simulation metadata explicitly unavailable. The post-hoc tables contain only fields consumed by public analyses. Prepared paired structures are provided under stable lowercase names.

`metadata/asset-rights.tsv` is the redistribution gate. Raw service packages, full molecular-dynamics artifacts, checkpoints, and licensed binary formats are not Git payloads; their availability is represented only by stable labels in `external_archive_manifest.tsv`. `metadata/manual_decisions.tsv` records reviewed corrections, annotations, parameter freezing, and representative-case selections.

`manifest.tsv` records every local file in this directory except itself and `SHA256SUMS`; the checksum file then covers every other file, including both manifests. Empty numeric fields mean unavailable source values. Tables use UTF-8, LF endings, lowercase booleans, and deterministic stable sorting.
