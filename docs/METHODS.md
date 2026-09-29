# AnnRE 0.1.0a1

A working, candidate-focused tool for testing gene annotation against long-read RNA evidence. It flags reference joins worth inspecting, discovers repeated connections between separate reference genes, and discovers internal exon splice contexts missing from a gene's reference models. Alternate GTFs supply exact-coordinate comparisons and review tracks.

The output is a **review queue with auditable evidence**, not a repaired annotation. This is an early research prototype: no trained classifier, calibrated error probability, or whole-genome performance claim. The BSF configuration demonstrates the tool on selected loci; it is not an independent benchmark.

## Install and run

Python 3.10+ and pysam are required. No pandas, plotting stack, browser service, or external annotation programs are needed.

```bash
python -m pip install /path/to/annre
annre prepare --config config.json --out run
annre extract --run run --library 0
# Repeat extract for every zero-based row in the library manifest.
annre report --run run --max-plots 30
```

Alternatively, from the source directory use `python -m annre`. `prepare` freezes the inputs, defaults, selected models and code checksums. A library can run independently in an array; `report` refuses missing or mismatched evidence. To retry a failed library, remove only that library's incomplete `run/evidence/NNNN` directory. Completed libraries are not silently overwritten. Use a new run directory after changing code or configuration.

Example configuration:

```json
{
  "genome": "genome.fa",
  "reference": "reference.gtf",
  "samples": "samples.tsv",
  "genes": ["gene_A", "gene_B"],
  "alternates": {"Tiberius": "predictions.gtf", "Bambu": "bambu.gtf"},
  "observation_mode": "pacbio-zmw",
  "strand_mode": "alignment"
}
```

Paths in the configuration are relative to its directory. BAM paths in the sample sheet are relative to the sample sheet. FASTA requires `.fai`; coordinate-sorted BAMs require indexes and matching reference contig names/lengths. Dictionary validation detects incompatible assemblies but does **not** prove sequence identity. Reference GTF input requires stranded exon records with gene/transcript IDs; unstranded alternate models are excluded explicitly and their counts are reported in the manifest and HTML; CDS-only GTF and GFF3 require explicit conversion upstream. Discontinuous CDS alone does not define transcript UTRs.

The required sample columns are `library_id`, `sample_id`, `bam`. Optional `tissue`, `stage`, `instar` propagate to support tables. One row is one library; repeat a sample ID for its top-up libraries. Distinct sample IDs are not automatically independent biological replicates.

`genes` selects reference genes to audit. Their spans plus `context_bp` (default 5,000) define fetch windows. Neighboring gene models provide context, so discovered connections can include neighbors. An empty/omitted gene list selects the full annotation, but this implementation has **not** been optimized or benchmarked for whole-genome throughput. Start with a candidate panel. The supplied `examples/config.json` is a portable template. Supply your own genome, annotation, indexed alignments and metadata.

For SLURM, submit one array index per library, with the report dependent on successful completion of the whole array. See `examples/slurm.sh` for a generic array template. A 72-library array uses `--array=0-71`, without a concurrency throttle.

## What the evidence means

All TSV/audit intervals use **0-based, half-open coordinates**. SVG coordinate labels use 1-based inclusive display. Exons are stored in increasing genomic order on either strand; splice motifs are oriented to the transcript.

Default filters require primary mapped alignments, MAPQ ≥20 and not 255, no supplementary/secondary/duplicate/QC-fail flags or SA tag, clipping ≤10%, and an NM tag with edit fraction ≤5%. Gaps are CIGAR `N`; deletions are not introns. Exact observed junctions require 20 aligned reference bases immediately on each side. Filters are recorded in each library manifest. These defaults fit the present high-accuracy BSF data; they need validation before use on another protocol.

In `pacbio-zmw` mode, the counting key is sample ID plus movie/ZMW; multiple Kinnex segments from that observation do not become independent supporting reads for the same event. This is conservative evidence counting, **not** a count of original RNA molecules. `read-name` mode is available for other protocols, with its own duplicate-ID assumptions. Raw supporting alignment-segment counts and read IDs remain available. Conflicting strands within an observation are not separately adjudicated.

`alignment` strand mode assumes reads have been oriented to their transcripts, as in the BSF processed data. For unoriented minimap2 splice alignments, `ts` mode combines the alignment flag with the `ts` tag and rejects reads without a usable tag. See the [minimap2 cookbook](https://github.com/lh3/minimap2/blob/master/cookbook.md) for the tag convention. This strand mode has not been benchmarked on an ONT cohort.

### Reference junction audit

Each event is a unique annotated intron within a gene, with its transcript-specific adjacent exon pairs. A flank is expressed when an accepted alignment overlaps at least 50 bases of that exon (or its full length if shorter).

- `EXACT_JUNCTION_SUPPORTED`: ≥3 deduplicated reads from ≥2 sample IDs support the exact junction with strict anchors.
- `ALTERNATIVE_CONNECTION_SUPPORTED`: at least one distinct overlapping alternative junction has ≥3 reads from ≥2 samples on reads that cover both immediate reference flanks, or ≥50 reference-exonic bases in the gene context on each side of the reference intron. The latter detects alternative connections that skip an unsupported reference exon. Immediate-flank and gene-context bridges are reported separately. Different alternative coordinates are not pooled to pass this threshold. This does not establish a small boundary correction; large rewiring is included.
- `ALTERNATIVE_JUNCTION_SUPPORTED_CONNECTION_UNRESOLVED`: a distinct strictly anchored junction shares a reference donor or acceptor and has ≥3 reads from ≥2 samples, but repeated evidence connecting the gene context on both sides is missing. An alternative transcript start inside an intron can produce this pattern. Counts appear in `best_shared_boundary_reads/samples`; these reads do not validate the full gene connection.
- `CONNECTION_EVIDENCE_REVIEW`: some accepted exact or bridging evidence exists but does not satisfy the repeated-junction criteria. A bridge can be retained-intron or other structure.
- `UNSUPPORTED_HIGH_INFORMATION`: no accepted exact, immediate-flank bridging or gene-context bridging evidence, ≥10 reads on each flank in aggregate, and both flanks detected in ≥2 common sample IDs. This is a split-review candidate, not proof that the gene is wrong.
- `UNSUPPORTED_LOW_INFORMATION`: evidence is insufficient to make the stronger absence argument.

Reference-event `supporting_reads` counts **exact** junction support. Thus zero is expected for the focal unsupported split candidate even when both sides are strongly expressed. Separate columns show left/right, immediate-flank bridges, gene-context bridges, best gene-connecting alternative, and best shared-boundary alternative support. `alternative_connections.tsv` records exact alternative coordinates, counts and evidence scope. A coordinate may occur in both scopes; those rows must not be added together.

### Cross-gene junctions

A strictly anchored junction connects read segments overlapping exons of different same-strand reference genes. Each gene also contributes ≥50 aligned exonic bases to the read. Repeated support (≥3 reads/≥2 samples) creates a review candidate. Overlapping gene annotations can make ownership ambiguous; these are not automatic merges.

`PRESERVES_BOTH_MODELS_READTHROUGH_COMPATIBLE` additionally requires repeated reads containing a complete annotated splice chain from each gene, with terminal exon overlap. This is compatible with readthrough and preserves a useful counterexample to automatic merging. It does not establish complete transcript ends, intact coding sequences, or the biological mechanism.

### Novel internal exon contexts

A read internal exon is considered missing from a gene's reference models when its exact interval is absent there and its immediate flanking read segments overlap that gene's annotated exons. Both flanking introns must meet the anchor threshold. Boundary variants, retained intronic sequence and larger exon changes can all enter this broad candidate class; it is not a count of newly discovered functional exons. The candidate identity includes the exon and **both exact introns**. Support is never pooled across different splice contexts. An exon may already occur in another gene's annotation; “novel” here is relative to the assigned reference gene.

Only contexts with immediate same-gene flank assignment are discovered in this release. The earlier BSF exploratory queue also contained unassigned and weakly anchored events; this output will not reproduce that larger queue one-for-one.

GT-AG, GC-AG and AT-AC are treated as canonical. Repeated noncanonical events remain explicitly flagged. Repetition alone cannot rule out systematic alignment or library artifacts.

### Alternate annotations

`alternate_SOURCE` lists models with the exact candidate intron, or both consecutive introns for a novel exon context. This is structural agreement, **not RNA evidence**. A final assembly that retains its supplied reference can therefore match an unsupported reference junction. Tiberius coding-only predictions are not evidence about UTR completeness. Model IDs are source-specific; a different ID does not imply a different splice chain.

## Outputs

- `candidates.tsv`: reference tests and discovered context-specific events, counts, categories and alternate model matches.
- `sample_support.tsv`: counts by event, evidence type, sample and tissue/stage. Alternative-junction rows include their exact coordinate pair in `detail`. Omitted rows mean no passing support, not absence of expression or a differential-expression result.
- `alternative_connections.tsv`: exact alternative coordinates and motifs for reference tests.
- `index.html` and `plots/*.svg`: offline review report and scalable figures. These are annotation/read-chain schematics, not IGV screenshots. Top three observed chains are grouped by exact intron chain; their ends are representative examples. Counts are event-associated and a read can support several events.
- `evidence/NNNN/observations.jsonl.gz`: arrays of `[event_id, evidence_type, observation_id, read_id, detail, chain_id]`. `detail` holds alternate coordinates or preserved reference model IDs where applicable. Library/sample identity comes from the manifest. Labels can produce multiple audit records for one read.
- `evidence/NNNN/chains.json`, `events.json`, `done.json`: chain geometry, event definitions, quality-filter counts and input/output provenance.
- `prepared.json`, `manifest.json`: frozen configuration, reference/alternate checksums, BAM size/mtime, FASTA-index checksum, software and code hashes. BAMs and the whole genome are not rehashed by this tool.

The counts are **events**, not distinct misannotated genes. Nearby novel exon boundaries can yield several events; no edit should be inferred by summing categories.

## Validation and contribution boundary

Run `python -m unittest discover -s tests -v` from the repository root. Tests exercise true and unsupported joins, low-information negatives, complete-model readthrough, alternative starts that do not validate gene connections, novel splice contexts, top-up deduplication, minus-strand motif orientation, GTF quoted semicolons, CIGAR deletions, quality rejection, incomplete arrays and evidence tampering.

The prospective contribution is an explicit layer for deciding **which existing annotation structures deserve review and why**, using long reads as evidence and alternate annotations as context. Automatic split/merge writing, protein homology/CDS checks, whole-genome indexing, and independently labelled precision/recall evaluation remain future work. Existing BSF adjudication is a useful demonstration set, not a held-out validation set.

No open-source license has been selected for this preview. See [license status](../LICENSE_STATUS.md).
