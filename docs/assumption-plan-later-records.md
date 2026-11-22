# Later-records test: counts

Block T1 of `docs/assumption-plan.md` holds the rule for the later-records test. This page holds the counts behind it. They were made on 21 Nov 2026 with no fit, and they are the only numbers taken from the later records before the scoring.

## File

- `later-records.jsonl` holds the 252,148 BV-BRC lab records inserted at or after 2026-11-05T19:10:00, on 20,108 genomes. Its sha256 is `2efc5328ce4267a7963b93578899545758035d4b6bf004034ac41596d3f6a7dc`.
- The genomes' metadata rows are beside it in `later-genomes.jsonl`, sha256 `592a806aef1ca006c1bad54aaa8419ea7ebecbaac8cfa5907bd5d672972fe6c6`.
- Both are cut from the download of 21 Nov 2026, 22:05:58Z to 22:10:40Z: 1,537,259 records. The 1,285,111 records before the cut are the training table.
- The files were closed on 21 Nov 2026 at 23:18:12Z. A gzip copy of the records has sha256 `e53ce1f7cc47bdba592f9e6dd386b81065bef15e655a88c15e7a111733c40b57`.

## Unseen

Under today's label rule the later records add 155,448 R or S labels. They partition into 101,593 in 166 of the 274 datasets and 53,855 outside the 274. The table's four rows partition the same labels by genome.

| Genome | Later labels | In the 274 | Outside |
|---|---:|---:|---:|
| In its own dataset's report | 0 | 0 | 0 |
| In another report of the species | 0 | 0 | 0 |
| In no report, with a lab record before the cut | 102 | 86 | 16 |
| No lab record before the cut | 155,346 | 101,507 | 53,839 |

- The 102 labels are each for a drug the genome had no label for before the cut. They stay in.
- 3 S labels in S. enterica ceftriaxone share a BioSample with a genome in that dataset's own report, where it has the same label. They stay out.
- The accession check reaches 19,343 of the release's 23,757 genomes. A copy of an isolate with no accession is not found this way.

## Scorable

Each count is of the 274 datasets with at least that many later labels in each class.

| Label rule | R | S | At 1 | At 10 | At 25 |
|---|---:|---:|---:|---:|---:|
| Today's | 43,599 | 57,994 | 142 | 116 | 102 |
| Mixed pairs set aside | 43,554 | 57,990 | 142 | 116 | 102 |
| Most resistant call | 43,671 | 57,990 | 142 | 116 | 102 |
| Third class | 43,554 | 57,990 | 142 | 116 | 102 |

Outside the 274, 630 species and drug pairs hold a later label. 74 of them are curated datasets on today's table, and 56 reach 25 in each class. No release model exists for them.

## Old labels

The training table holds 612,989 labelled pairs. 9 gain a later record, and none changes its mix of calls. Under each of the four rules no old label changes or is lost.

## Genomes

The labels that stay in are on 12,617 genomes with 54.85 Gbp of assembly. 10,813 of those genomes, 48.11 Gbp, hold a label in the 274. The largest species are S. enterica (3,894), K. pneumoniae (2,502) and M. tuberculosis (1,192). No sequence is fetched before the scoring.

## Models

The release publishes its SCM and CART models, and its lasso when every feature is listed. Those are scored as they stand. It publishes no forest or boosted-tree model (`docs/prediction-contract.md`). For those two learners the default arm's refit saves its final models, and it reads no later record. A refit is a new model, so its rows are named as the refit and never as the release's model. Each arm saves its final models the same way.
