# DSA4262 Long-read RNA-seq Workflow

This repository contains the Nextflow workflow used for Task 5 of the DSA4262 project.

The workflow was adapted from the long-read RNA-seq pipeline introduced in the course workshop.

## Workflow

The pipeline performs four main steps:

1. Read alignment using Minimap2
2. SAM-to-BAM conversion using samtools
3. Transcript discovery and quantification using Bambu
4. Quality control using NanoPlot and samtools stats

The workflow supports both direct RNA and cDNA samples and uses different Minimap2 options for the two read types.

## Task 5 scenarios

Two Bambu scenarios were tested:

- Scenario 1: Bambu with genome annotations
- Scenario 2: Bambu without genome annotations using Nextflow `-resume`

The input paths in `task5_pipeline.nf` correspond to the directory structure used for the assignment.
