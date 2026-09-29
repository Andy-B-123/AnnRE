# Development and attribution

**Generated with OpenAI Codex and OpenAI models; reviewed and guided by Andreas Bachler.**

Codex generated and revised substantial parts of the implementation, tests, documentation, analysis workflows and presentation assets. Andreas Bachler supplied the biological questions, datasets, requirements, scientific interpretation, review and iterative guidance. This is an early research prototype: that review does not imply that every generated line has undergone independent code review or that candidate classifications have been independently validated.

The project maintainer is responsible for release decisions, interpretation and suitability for use. Automated tests and selected BSF regression analyses are described separately from biological validation. Codex is acknowledged as a development tool; it is not listed as a human scientific author. This project is not endorsed by OpenAI.

Model versions varied during development and have not been exhaustively verified from session logs. No exact model identifier is asserted here. OpenAI's Codex documentation: https://developers.openai.com/codex/

## Scientific software and data

AnnRE uses pysam to read indexed alignments and FASTA files. Reference annotations, genome assemblies, alignment workflows and alternate annotations require their own provenance and citations. The BSF demonstration used NCBI RefSeq, Tiberius predictions, Bambu output and an existing ANNEXA trial. AnnRE does not bundle or replace those tools.

The published comparative annotation work is separate from the BSF prototype: https://doi.org/10.1186/s12864-025-11765-w . Its 605 cases across 30 genomes are not results of this software release.

Add agreed scientific contributors, data-generation credits and third-party software citations before a formal research release. The current citation file names the project lead only and should be updated as authorship is agreed.
