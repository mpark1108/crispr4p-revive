# PomBase annotation snapshots

CRISPR4P ships fixed PomBase datasets so name lookup and cut-site annotation
are reproducible and work from a fresh clone without network access.

## Gene names and synonyms

- File: `gene_IDs_names.tsv`
- Purpose: current systematic IDs, primary gene names, and searchable synonyms
- PomBase source directory:
  <https://www.pombase.org/data/names_and_identifiers/>
- Embedded Chado database date: 2026-08-24 20:37
- SHA-256:
  `4c688312ebb5ab80356cf60c0b85ad9f409d049e5688c0915e477de60cde4899`

The smaller names-only table is packaged instead of
`gene_IDs_names_products.tsv` because both contain the same name and synonym
fields and CRISPR4P does not currently display gene-product descriptions.

## Genome feature annotations

- File: `Schizosaccharomyces_pombe_all_chromosomes.gff3`
- Purpose: gene, transcript, CDS, intron, UTR, and non-coding-exon coordinates
- PomBase source directory:
  <https://www.pombase.org/data/genome_sequence_and_features/gff3/>
- Local source-file timestamp: 2026-07-30
- SHA-256:
  `88e4a26c16762c6d97e7f6a1600cd7dac193bdaeff4753cc1dd2123a79f6a025`

## Gene viability summary

- File: `gene_viability.tsv`
- Purpose: PomBase null/deletion viability summary by systematic gene ID
- PomBase documentation:
  <https://www.pombase.org/downloads/phenotype-annotations>
- Local source-file timestamp: 2026-08-07
- SHA-256:
  `e9399024327be0a2a6618c8fda1dfaef6ef72c493eb1050c7a6d86d41a3d3d09`

The viability file contains the raw states `viable`, `inviable`,
`condition-dependent`, and `unknown`.

The GFF3 and viability files do not embed an exact PomBase release identifier.
Their timestamps are those of the reviewed local copies. The gene-names table
does embed its Chado database date. The hashes are the authoritative
identifiers for every snapshot used by this project.
