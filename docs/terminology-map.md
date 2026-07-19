# Terminology map

| Manuscript display name | Code/config name | Principal data column | Meaning in this repository |
|---|---|---|---|
| OIPS-P_static | `OIPS-P_static` | `OIPS-P_static` | Frozen static candidate-priority score. |
| OIPS-P | contextual/legacy label | `OIPS-P_static` when explicitly static | Do not use without a qualifier where dynamic or confidence layers could be confused. |
| O_rel | `O_rel_formal` | `O_rel_formal` | Static oligomeric-relevance module calculated from the prepared assembly proxy. |
| O_rel_formal | `O_rel_formal` | `O_rel_formal` | Code/data serialization of manuscript O_rel. |
| Q_evidence | `Q_evidence` | `Q_evidence` | Static evidence-completeness/quality module; not post-hoc validation confidence. |
| A_rob | not implemented as a released independent module | none | Historical or manuscript term requiring author confirmation; it must not be silently equated with Q_evidence. |
| OIPS-T | not implemented | none | Pocket-type interpretation layer outside the released static score. |
| OIPS-C | not implemented | none | Evidence-confidence layer outside the released static score. |
| reference-associated | reference endpoint | reference mapping and target metrics | First static-ranked candidate meeting the frozen reference rule. |
| first-supported | first-supported endpoint | target metrics | First candidate supported by the released post-static evidence rules; not independent validation. |
| assembly-supported | evidence label/context | automated evidence labels | Candidate with assembly-relevant support; not synonymous with reference-associated. |

Existing public column names are retained to preserve frozen output identity.
Potential conceptual overlaps are documented rather than bulk-renamed.
