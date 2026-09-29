# BSF evidence summary

The current refined shortlist contains **105 genes with 116 suspect reference joins**, about **0.75% of the 13,912 intron-containing reference genes** assessed. This is a review burden, not an estimated true annotation-error rate.

The project-specific genome-wide census examined 77,385 unique reference introns in 72 long-read libraries grouped into 68 sample IDs. The reference has 15,712 genes; the 1,800 without reference introns were outside this intron screen. The portable AnnRE alpha was demonstrated on a ten-gene panel; it has not been benchmarked genome-wide.

| Result | Count |
|---|---:|
| Refined review candidates | 105 genes / 116 joins |
| Subset with nearby alignment-end evidence | 55 genes / 56 joins |
| Reference transcript models containing the 116 joins | 160 |
| Suspect joins retained by Bambu (72 libraries) | 116 / 116 |
| Suspect joins retained by trial Bambu (12 libraries) | 116 / 116 |
| Suspect joins retained by ANNEXA final output (12 libraries) | 116 / 116 |
| Separate alternate gene models on both sides, Bambu 72-library output | 1 / 116 events |
| Separate alternate gene models on both sides, ANNEXA final output | 4 / 116 events |
| Tiberius separate models on both sides without a model bridging both anchor regions | 61 / 116 events |

The Bambu and ANNEXA outputs also retain all 160 corresponding reference transcript models. Adding alternative models and removing an unsupported reference structure are different operations. The comparison measures **candidate retention**, not overall annotation quality or accuracy. AnnRE itself currently flags candidates; it does not remove these joins automatically.

Of the 116 events, 67 joins in 62 genes meet the project's expression-information rule in the ANNEXA 12-library subset; all 67 remain in its annotation output. Selection used all 72 libraries, so this is not an independent held-out benchmark. Bambu and ANNEXA were not evaluated with equivalent full-cohort inputs.

Bambu assigns nonzero abundance to 149 of the 160 retained reference models, but its full-length count is zero for every one. Nonzero assigned abundance does not establish that a read supports the complete transcript structure. Tiberius predictions are complementary coding-model evidence, not biological ground truth.

## Why the count changed from 108 to 105

The initial 50-base anchor screen shortlisted 108 genes / 119 joins. Follow-up inspection recovered shorter-anchor exact junction support in three genes: LOC119661631 (70 observations meeting a 20-base anchor rule across five sample IDs), LOC119655972 (two across two samples), and LOC119651052 (one in one sample). These were removed from the no-observed-exact-support shortlist. All 72 libraries were included in the follow-up anchor audit.

The remaining candidates still need adjudication. Low expression, alignment limitations, partial transcripts, readthrough and overlapping transcription can complicate interpretation. Alternative bridging structures have not been exhaustively reassessed with every possible short-anchor rule.

## Suggested conference wording

> A genome-wide long-read evidence audit shortlisted 105 BSF genes (116 reference joins; approximately 0.75% of intron-containing genes) for annotation review. These joins persisted in existing Bambu and ANNEXA outputs, illustrating a complementary role for evidence-based reference-model QC. The candidates require biological and structural adjudication; this is not a measured accuracy improvement.

The archived project evidence includes per-junction comparisons, read audits, quantification and checksums. This source repository contains a portable prototype and a summary, not the underlying BAMs or the entire project-specific census. Counts alone cannot reproduce that census. Independent precision/recall or validated correction yield remains future work.
