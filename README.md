# RNA-seq Diagnostic Plots: A Practical Guide

A reference guide for quality control in bulk RNA-seq analysis, organized in workflow order from raw sequencing quality through outlier detection. Each plot is annotated with what good looks like, what bad looks like, and what action to take.

Plots 01-02 use simulated data to illustrate patterns that require raw FASTQ or BAM files in practice. Plots 03-06 use real data from GSE161731 downloaded directly from GEO.

## Biological Question

How do you know if your RNA-seq data is high quality before running differential expression analysis? This repo provides a six-plot QC workflow that catches the most common problems: failed sequencing runs, contamination, outlier samples, and biotype composition issues before they affect downstream results.

## Disclaimer

This is an independent secondary analysis of publicly available data from GEO (GSE161731) for educational and portfolio purposes. The original study was designed and conducted by McClain et al. (2021).

## Dataset

**GSE161731**: Blood RNA-seq from patients with respiratory infections and healthy controls. McClain et al. *Nature Communications* 12, 1079 (2021). https://doi.org/10.1038/s41467-021-21289-y

Data downloads automatically from GEO on first run; no manual download required.

## Plots Covered

### Stage 1: Pre-alignment QC (raw FASTQ)

> These plots require raw FASTQ files in practice. Simulated data illustrates typical patterns. Generate with FastQC and aggregate with MultiQC.

| Plot | Data | What it shows | Tool (in practice) |
|---|---|---|---|
| 01. Per-base sequence quality | Simulated | Phred quality scores across 150bp read positions; good vs 3' degradation | FastQC + MultiQC |
| 02. Genomic feature assignment | Simulated | Exon/intron/intergenic/rRNA fractions; good vs four failure modes | featureCounts (Subread) |

### Stage 2: Count-level QC (after alignment and counting)

| Plot | Data | What it shows |
|---|---|---|
| 03. Library size distribution | Real | Total mapped reads per sample, used to flag low-depth samples |
| 04. Gene detection per sample | Real | Number of genes detected per sample, to cross-reference with library size |
| 05. RNA biotype distribution | Real | Proportion of reads mapping to protein-coding, ncRNA, rRNA, etc. |

### Stage 3: Sample-level QC (outlier detection)

| Plot | Data | What it shows |
|---|---|---|
| 06. PCA outlier detection | Real | Identifies extreme samples more than 3 SD from the mean on PC1 or PC2 |

## Key Interpretations

### Per-base sequence quality (Plot 01)
Good: Phred scores consistently above Q30 with a slight acceptable decline at the 3' end. Bad: sharp drop below Q20 after position 100: trim with Trimmomatic (`SLIDINGWINDOW:4:20`) or Cutadapt (`--quality-cutoff 20`) before alignment.

### Genomic feature assignment (Plot 02)
Good: more than 60% exonic reads. High intronic reads (more than 40%) suggest genomic DNA contamination. High intergenic reads (more than 30%) suggest severely degraded RNA or wrong genome annotation. High rRNA (more than 20%) indicates failed rRNA depletion.

### Library size distribution (Plot 03)
Good: consistent depth across samples, median 20-80M reads for blood RNA-seq. Flag samples below 20% of the median (dotted red line). Always cross-reference with gene detection: low library size AND low gene detection is a strong indicator of a failed sample.
Note: One CoV_other sample shows unusually high library size (~130M reads, approximately 2.5x the median). This sample passes gene detection and PCA QC checks and was retained: high library size alone is not grounds for removal without additional evidence of quality failure.
### Gene detection per sample (Plot 04)
Good: 12,000-25,000 genes detected per sample for blood RNA-seq, consistent within cohorts. Samples below 5,000 detected genes should be investigated and likely removed.

### RNA biotype distribution (Plot 05)
Good: more than 70% protein-coding (dashed line), consistent across all samples. High rRNA suggests failed rRNA depletion. A sample with an unusual biotype profile relative to its cohort mates may be mislabeled or from a failed library prep.

### PCA outlier detection (Plot 06)
Good: samples cluster loosely by biological group, no single sample isolated far from the rest. Samples more than 3 SD from the mean on PC1 or PC2 (red circles) are candidates for removal; always cross-reference with library size and gene detection before deciding. In this dataset three outliers were identified and removed before downstream analysis.

## What This Dataset Shows

GSE161731 is a high-quality clinical RNA-seq dataset that passes all QC checks:

- Library sizes are consistent (median 51M reads, no samples below threshold)
- Gene detection is consistent across cohorts (15,000-25,000 genes per sample)
- RNA biotype composition is dominated by protein-coding genes (more than 95% across all samples)
- Three outlier samples were correctly identified by PCA and removed before analysis

This makes it a useful reference for what good RNA-seq QC looks like in practice.

## Tools

**R:** edgeR, tidyverse, ggplot2, ggrepel, org.Hs.eg.db, AnnotationDbi, RColorBrewer

**External (for pre-alignment QC, not in script):** FastQC, MultiQC, featureCounts (Subread)

## Project Structure

```
scripts/
└── 01_rnaseq_diagnostic_plots.R    # All 6 plots, self-contained

figures/
├── 01_per_base_quality_simulated.png
├── 02_genomic_feature_assignment_simulated.png
├── 03_library_size.png
├── 04_gene_detection.png
├── 05_biotype_distribution.png
└── 06_pca_outlier_detection.png

data/                               # Downloaded from GEO on first run (not tracked)
```

## How to Reproduce

1. Clone this repository
2. Open `scripts/01_rnaseq_diagnostic_plots.R` in RStudio
3. Run from top to bottom, and the data downloads automatically on first run
4. All 6 figures saved to `figures/`

## Limitations

- Plots 01 and 02 use simulated data. In practice these require raw FASTQ files (FastQC/MultiQC) and alignment BAM files (featureCounts/RSeQC).
- Biotype annotations use org.Hs.eg.db which maps Ensembl IDs to broad gene type categories. More granular biotype annotations are available from Ensembl BioMart.
- QC thresholds used here (20% of median library size, 5,000 gene minimum) are appropriate for blood RNA-seq and may need adjustment for other tissue types.
- All real-data plots use post-QC samples (n=138 after outlier removal and deduplication). Plots 03 and 04 intentionally use pre-QC data to demonstrate what outliers look like.

## References

McClain, M.T., Constantine, F.J., Henao, R. et al. Dysregulated transcriptional responses to SARS-CoV-2 in the periphery. *Nature Communications* 12, 1079 (2021). https://doi.org/10.1038/s41467-021-21289-y

Andrews, S. FastQC: A quality control tool for high throughput sequence data. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Liao, Y., Smyth, G.K. & Shi, W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics* 30, 923-930 (2014). https://doi.org/10.1093/bioinformatics/btt656

## License

MIT License. See [LICENSE](LICENSE) for details.
