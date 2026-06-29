# ==============================================================================
# RNA-seq Diagnostic Plots: A Practical Guide
# Dataset: GSE161731 — Blood Transcriptomics in Respiratory Infections
#
# This script generates key diagnostic plots used in bulk RNA-seq analysis,
# organized in workflow order from raw sequencing quality through outlier
# detection. For each plot the code documents:
#   - What the plot shows
#   - What a good plot looks like
#   - What a bad plot looks like and what it suggests
#   - What action to take
#
# Plots 01-02 use simulated data (require raw FASTQ/BAM files in practice).
# Plots 03-06 use real data from GSE161731 downloaded directly from GEO.
#
# Plots covered:
#   01. Per-base sequence quality    (simulated — FastQC)
#   02. Genomic feature assignment   (simulated — featureCounts)
#   03. Library size distribution    (real data)
#   04. Gene detection per sample    (real data)
#   05. RNA biotype distribution     (real data)
#   06. PCA — outlier detection      (real data)
# ==============================================================================


# ==============================================================================
# 1. SETUP
# ==============================================================================

library(edgeR)
library(tidyverse)
library(ggrepel)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(RColorBrewer)

setwd("~/Desktop/rnaseq-diagnostic-plots")
dir.create("figures", showWarnings = FALSE)

palette_cohort <- c(
  healthy   = "#2ecc71",
  COVID19   = "#e74c3c",
  Influenza = "#3498db",
  Bacterial = "#e67e22",
  CoV_other = "#9b59b6"
)


# ==============================================================================
# 2. DOWNLOAD DATA FROM GEO
# ==============================================================================
#
# GSE161731: Blood RNA-seq from COVID-19, Influenza, Bacterial, CoV_other,
# and Healthy controls. McClain et al. Nature Communications 2021.
#
# NOTE: GSE161731 stores counts in supplementary files, not the series matrix.
# The key file maps RNA sample IDs to metadata (cohort, age, gender, batch).
install.packages("data.table")
cat("Downloading GSE161731 from GEO...\n")
dir.create("data", showWarnings = FALSE)

count_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161731/suppl/GSE161731_counts.csv.gz"
key_url   <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161731/suppl/GSE161731_counts_key.csv.gz"

if (!file.exists("data/GSE161731_counts.csv.gz"))
  download.file(count_url, "data/GSE161731_counts.csv.gz", mode = "wb")
if (!file.exists("data/GSE161731_counts_key.csv.gz"))
  download.file(key_url, "data/GSE161731_counts_key.csv.gz", mode = "wb")

counts_raw <- read.csv("data/GSE161731_counts.csv.gz", row.names = 1)
key        <- as.data.frame(data.table::fread("data/GSE161731_counts_key.csv.gz"),
                            stringsAsFactors = FALSE)

# Fix column name mismatch — R converts dashes to dots and adds X prefix
key$rna_id_fixed <- make.names(key$rna_id)
rownames(key)    <- key$rna_id_fixed

overlap <- intersect(colnames(counts_raw), key$rna_id_fixed)
counts  <- as.matrix(counts_raw[, overlap])
meta        <- as.data.frame(key[overlap, ], stringsAsFactors = FALSE)
meta$age    <- as.numeric(meta$age)
meta$batch  <- as.numeric(meta$batch)

cat("Counts:", nrow(counts), "genes x", ncol(counts), "samples\n")
cat("Metadata:", nrow(meta), "samples\n")

# Remove all-zero genes
counts_filtered <- counts[rowSums(counts) > 0, ]
cat("After removing zero-count genes:", nrow(counts_filtered), "genes\n")

# Clean cohort names
meta$cohort <- recode(meta$cohort,
                      "COVID-19"  = "COVID19",
                      "CoV other" = "CoV_other")
meta <- meta[colnames(counts_filtered), ]


# ==============================================================================
# 3. QC — REMOVE OUTLIERS AND APPLY GENE FILTER
# ==============================================================================
# Applied before plots 05-06 which require cleaner data.
# Plots 03-04 use counts_filtered (pre-QC) to show what outliers look like.

# Remove extreme outliers identified by PCA (> 3 SD on PC1 or PC2)
dge_raw  <- DGEList(counts = counts_filtered)
log_cpm  <- cpm(dge_raw, log = TRUE, prior.count = 1)
pca_raw  <- prcomp(t(log_cpm), scale. = FALSE)

pca_df <- data.frame(
  PC1    = pca_raw$x[, 1],
  PC2    = pca_raw$x[, 2],
  sample = colnames(counts_filtered),
  cohort = meta$cohort
)

pc1_mean <- mean(pca_df$PC1); pc1_sd <- sd(pca_df$PC1)
pc2_mean <- mean(pca_df$PC2); pc2_sd <- sd(pca_df$PC2)
pca_df$outlier <- abs(pca_df$PC1 - pc1_mean) > 3 * pc1_sd |
  abs(pca_df$PC2 - pc2_mean) > 3 * pc2_sd

outlier_samples <- pca_df$sample[pca_df$outlier]
cat("Outliers removed:", length(outlier_samples), "\n")

counts_qc <- counts_filtered[, !colnames(counts_filtered) %in% outlier_samples]
meta_qc   <- as.data.frame(meta[colnames(counts_qc), ],
                           stringsAsFactors = FALSE)
# Remove missing gender and age
meta_qc <- meta_qc %>%
  filter(gender != "" & !is.na(gender) & !is.na(age))
meta_qc <- as.data.frame(meta_qc)
counts_qc <- counts_qc[, rownames(meta_qc)]

# Keep one sample per subject (earliest timepoint)
time_rank <- c("early" = 1, "middle" = 2, "late" = 3)
meta_qc$time_rank <- ifelse(
  is.na(meta_qc$time_since_onset) | meta_qc$time_since_onset == "",
  0, time_rank[meta_qc$time_since_onset]
)
meta_qc <- meta_qc %>%
  group_by(subject_id) %>%
  arrange(time_rank) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  as.data.frame()
meta_qc <- as.data.frame(meta_qc)
rownames(meta_qc) <- meta_qc$rna_id_fixed
counts_qc <- counts_qc[, rownames(meta_qc)]

# Gene filtering — keep genes detected in at least 10% of samples
counts_qc <- counts_qc[rowSums(counts_qc >= 10) >= 0.1 * ncol(counts_qc), ]

meta_qc$cohort <- factor(meta_qc$cohort,
                         levels = c("healthy","COVID19","Influenza",
                                    "Bacterial","CoV_other"))
cat("After QC:", nrow(meta_qc), "samples,", nrow(counts_qc), "genes\n")


# ==============================================================================
# PLOT 01: PER-BASE SEQUENCE QUALITY (SIMULATED)
# ==============================================================================
#
# WHAT IT SHOWS:
#   Phred quality scores at each position along the sequencing read.
#   Q = -10 * log10(error probability). Q30 = 1 error per 1,000 bases.
#   Generated by FastQC; aggregated across samples with MultiQC.
#
#   Simulated here because this plot requires raw FASTQ files.
#   In practice run FastQC before alignment and inspect before proceeding.
#
# WHAT GOOD LOOKS LIKE:
#   - Scores above Q30 across all positions
#   - Slight decline toward the 3' end is normal if it stays above Q28
#
# WHAT BAD LOOKS LIKE:
#   - Sharp quality drop after position 100 (3' end degradation)
#     → trim with Trimmomatic (SLIDINGWINDOW:4:20) or Cutadapt
#   - Low quality throughout (< Q20) → sequencing run may have failed
#   - Sudden drop at one position → possible instrument cycle failure
#
# ACTION:
#   Trim low-quality 3' ends. Re-run FastQC after trimming to confirm.

set.seed(42)
positions   <- 1:150

good_mean <- 36 - (positions / 150) * 4
bad_mean  <- ifelse(positions <= 100,
                    35 - (positions / 100) * 5,
                    30 - ((positions - 100) / 50) * 18)

qual_df <- bind_rows(
  data.frame(position = positions,
             mean  = good_mean,
             q25   = good_mean - 2,
             q75   = good_mean + 1.5,
             quality = "Good quality"),
  data.frame(position = positions,
             mean  = bad_mean,
             q25   = bad_mean - 3,
             q75   = bad_mean + 2,
             quality = "Poor quality (3' degradation)")
)

p01 <- ggplot(qual_df, aes(x = position)) +
  annotate("rect", xmin=0, xmax=151, ymin=28, ymax=40,
           fill="#2ecc71", alpha=0.10) +
  annotate("rect", xmin=0, xmax=151, ymin=20, ymax=28,
           fill="#f39c12", alpha=0.10) +
  annotate("rect", xmin=0, xmax=151, ymin=0,  ymax=20,
           fill="#e74c3c", alpha=0.10) +
  annotate("text", x=148, y=38, label="Good (Q≥28)",
           color="#27ae60", size=3, hjust=1, fontface="bold") +
  annotate("text", x=148, y=24, label="Acceptable (Q20-28)",
           color="#e67e22", size=3, hjust=1, fontface="bold") +
  annotate("text", x=148, y=10, label="Poor (Q<20)",
           color="#c0392b", size=3, hjust=1, fontface="bold") +
  geom_ribbon(aes(ymin=q25, ymax=q75, fill=quality), alpha=0.25) +
  geom_line(aes(y=mean, color=quality), linewidth=1) +
  geom_hline(yintercept=c(20,30), linetype="dashed",
             color="grey50", linewidth=0.5) +
  scale_color_manual(values=c("Good quality"="#27ae60",
                              "Poor quality (3' degradation)"="#c0392b"),
                     name="") +
  scale_fill_manual(values=c("Good quality"="#27ae60",
                             "Poor quality (3' degradation)"="#c0392b"),
                    name="") +
  scale_y_continuous(limits=c(0,42), breaks=seq(0,40,10)) +
  scale_x_continuous(breaks=seq(0,150,25)) +
  facet_wrap(~quality, ncol=2) +
  labs(
    title    = "01. Per-Base Sequence Quality (Simulated)",
    subtitle = paste0("Simulated 150bp reads | Line = mean Phred score | Ribbon = IQR\n",
                      "In practice: generate with FastQC, aggregate with MultiQC"),
    x = "Position in Read (bp)",
    y = "Phred Quality Score"
  ) +
  theme_minimal(base_size=12) +
  theme(plot.title      = element_text(face="bold"),
        plot.subtitle   = element_text(size=9, color="grey40"),
        strip.text      = element_text(face="bold", size=11),
        legend.position = "none",
        panel.grid.minor = element_blank())

ggsave("figures/01_per_base_quality_simulated.png", p01,
       width=12, height=5, dpi=300)
cat("Plot 01 saved: per-base sequence quality\n")


# ==============================================================================
# PLOT 02: GENOMIC FEATURE ASSIGNMENT (SIMULATED)
# ==============================================================================
#
# WHAT IT SHOWS:
#   Proportion of reads assigned to exons, introns, intergenic regions,
#   rRNA, and unassigned. Generated by featureCounts (Subread) or RSeQC.
#
#   Simulated here because this plot requires aligned BAM files.
#
# WHAT GOOD LOOKS LIKE:
#   - More than 60% exonic reads
#   - Less than 20% intronic (some is normal — pre-mRNA)
#   - Less than 10% intergenic
#   - Low unassigned (less than 10%)
#
# WHAT BAD LOOKS LIKE:
#   - High intronic (more than 40%) → genomic DNA contamination
#   - High intergenic (more than 30%) → severely degraded RNA or wrong genome
#   - High rRNA (more than 20%) → rRNA depletion failed during library prep
#   - High unassigned (more than 30%) → wrong genome used for alignment
#
# ACTION:
#   High intronic: check RIN score, consider RNase H treatment.
#   High rRNA: repeat library prep with better rRNA depletion.
#   High unassigned: verify genome and annotation versions match.

scenarios <- list(
  "Good\n(typical blood RNA-seq)" = c(
    Exonic=72, Intronic=18, Intergenic=5, rRNA=2, Unassigned=3),
  "High Intronic\n(DNA contamination)" = c(
    Exonic=38, Intronic=48, Intergenic=7, rRNA=2, Unassigned=5),
  "High Intergenic\n(degraded RNA)" = c(
    Exonic=35, Intronic=20, Intergenic=38, rRNA=3, Unassigned=4),
  "High rRNA\n(failed rRNA depletion)" = c(
    Exonic=30, Intronic=10, Intergenic=5, rRNA=52, Unassigned=3)
)

feature_df <- map_dfr(names(scenarios), function(s) {
  vals <- scenarios[[s]]
  data.frame(scenario=s, feature=names(vals), pct=as.numeric(vals))
})

feature_df$feature  <- factor(feature_df$feature,
                              levels=c("Exonic","Intronic",
                                       "Intergenic","rRNA","Unassigned"))
feature_df$scenario <- factor(feature_df$scenario, levels=names(scenarios))

feature_colors <- c(Exonic="#2ecc71", Intronic="#3498db",
                    Intergenic="#e67e22", rRNA="#e74c3c",
                    Unassigned="#95a5a6")

p02 <- ggplot(feature_df, aes(x=scenario, y=pct, fill=feature)) +
  geom_bar(stat="identity", width=0.7) +
  geom_hline(yintercept=60, linetype="dashed",
             color="black", linewidth=0.6) +
  annotate("text", x=0.55, y=62,
           label="60% exonic threshold", size=3, hjust=0) +
  scale_fill_manual(values=feature_colors, name="Genomic Feature") +
  scale_y_continuous(breaks=seq(0,100,20), limits=c(0,105)) +
  geom_text(aes(label=paste0(pct,"%")),
            position=position_stack(vjust=0.5),
            size=3, color="white", fontface="bold") +
  labs(
    title    = "02. Genomic Feature Assignment (Simulated)",
    subtitle = paste0("Simulated scenarios showing common failure modes\n",
                      "In practice: generate with featureCounts (Subread) or RSeQC"),
    x = "", y = "Percentage of Reads (%)"
  ) +
  theme_minimal(base_size=12) +
  theme(plot.title       = element_text(face="bold"),
        plot.subtitle    = element_text(size=9, color="grey40"),
        axis.text.x      = element_text(size=10, face="bold"),
        panel.grid.minor = element_blank())

ggsave("figures/02_genomic_feature_assignment_simulated.png", p02,
       width=11, height=6, dpi=300)
cat("Plot 02 saved: genomic feature assignment\n")


# ==============================================================================
# PLOT 03: LIBRARY SIZE DISTRIBUTION
# ==============================================================================
#
# WHAT IT SHOWS:
#   Total mapped reads per sample (library size), ordered from lowest to
#   highest. Each bar is one sample, colored by cohort.
#
# WHAT GOOD LOOKS LIKE:
#   - Consistent depth across samples — median 20-80M reads for blood RNA-seq
#   - No samples more than 3-5x below the median
#
# WHAT BAD LOOKS LIKE:
#   - Samples below 20% of the median (dotted red line)
#     → failed sequencing, poor RNA quality, or sample swap
#   - Bimodal distribution → two different sequencing runs or batches
#   - Extreme outlier on the high end → possible sample duplication
#
# ACTION:
#   Flag samples below the red threshold. Cross-reference with gene detection
#   (Plot 04) — low library size AND low gene detection = strong QC failure.

lib_sizes <- colSums(counts_filtered) / 1e6

lib_df <- data.frame(
  sample   = names(lib_sizes),
  lib_size = lib_sizes,
  cohort   = meta$cohort
) %>% arrange(lib_size)
lib_df$sample <- factor(lib_df$sample, levels = lib_df$sample)

median_lib    <- median(lib_sizes)
low_threshold <- median_lib * 0.2

p03 <- ggplot(lib_df, aes(x=sample, y=lib_size, fill=cohort)) +
  geom_bar(stat="identity", width=0.8) +
  geom_hline(yintercept=median_lib, linetype="dashed",
             color="black", linewidth=0.6) +
  geom_hline(yintercept=low_threshold, linetype="dotted",
             color="red", linewidth=0.6) +
  scale_fill_manual(values=palette_cohort, name="Cohort") +
  annotate("text", x=5, y=median_lib+2,
           label=paste0("Median: ", round(median_lib,1), "M"),
           size=3.5, hjust=0) +
  annotate("text", x=5, y=low_threshold+2,
           label="Low threshold (20% of median)",
           size=3, hjust=0, color="red") +
  labs(
    title    = "03. Library Size Distribution",
    subtitle = paste0("GSE161731 | n=", ncol(counts_filtered),
                      " samples | Samples below red line flagged as low-depth"),
    x = "Sample (ordered by library size)",
    y = "Library Size (millions of reads)"
  ) +
  theme_minimal(base_size=12) +
  theme(axis.text.x     = element_blank(),
        axis.ticks.x    = element_blank(),
        plot.title      = element_text(face="bold"),
        legend.position = "right")

ggsave("figures/03_library_size.png", p03, width=12, height=5, dpi=300)
cat("Plot 03 saved: library size distribution\n")


# ==============================================================================
# PLOT 04: GENE DETECTION PER SAMPLE
# ==============================================================================
#
# WHAT IT SHOWS:
#   Number of genes with count > 0 per sample, grouped by cohort.
#
# WHAT GOOD LOOKS LIKE:
#   - 12,000-18,000 genes detected per sample for blood RNA-seq
#   - Consistent detection within cohorts (small variance)
#
# WHAT BAD LOOKS LIKE:
#   - Samples with fewer than 5,000 genes detected
#     → poor RNA quality, low sequencing depth, or failed extraction
#   - Unusually high gene detection → possible genomic DNA contamination
#   - High within-cohort variance → inconsistent sample handling
#
# ACTION:
#   Investigate samples below 5,000 detected genes. Always check alongside
#   library size — low library size AND low gene detection = failed sample.

gene_detect <- data.frame(
  sample  = colnames(counts_filtered),
  n_genes = colSums(counts_filtered > 0),
  cohort  = meta$cohort,
  lib_size = lib_sizes
)

p04 <- ggplot(gene_detect, aes(x=cohort, y=n_genes, color=cohort)) +
  geom_jitter(width=0.2, size=2, alpha=0.7) +
  geom_boxplot(alpha=0.3, outlier.shape=NA) +
  geom_hline(yintercept=5000, linetype="dashed",
             color="red", linewidth=0.5) +
  annotate("text", x=0.6, y=5300,
           label="5,000 gene minimum", size=3,
           color="red", hjust=0) +
  scale_color_manual(values=palette_cohort) +
  labs(
    title    = "04. Gene Detection per Sample",
    subtitle = "GSE161731 | Samples with fewer than 5,000 genes detected should be investigated",
    x = "", y = "Number of Genes Detected (count > 0)"
  ) +
  theme_minimal(base_size=12) +
  theme(plot.title      = element_text(face="bold"),
        legend.position = "none")

ggsave("figures/04_gene_detection.png", p04, width=8, height=5, dpi=300)
cat("Plot 04 saved: gene detection per sample\n")


# ==============================================================================
# PLOT 05: RNA BIOTYPE DISTRIBUTION
# ==============================================================================
#
# WHAT IT SHOWS:
#   Proportion of reads mapping to protein-coding genes, ncRNAs, pseudogenes,
#   rRNA, and other biotypes. Annotated using Ensembl gene biotypes from
#   org.Hs.eg.db. Samples ordered by cohort.
#
# WHAT GOOD LOOKS LIKE:
#   - More than 70% protein-coding for blood mRNA-seq with poly-A selection
#     or rRNA depletion (dashed line)
#   - Consistent proportions across samples within each cohort
#   - Low rRNA (less than 5%)
#
# WHAT BAD LOOKS LIKE:
#   - High rRNA (more than 20%) → rRNA depletion failed during library prep
#   - High pseudogene reads → possible genomic DNA contamination
#   - Sample with unusually different biotype profile vs cohort mates
#     → mislabeled sample or failed library prep
#   - Protein-coding consistently below 70% → check library prep protocol
#
# ACTION:
#   Flag samples with unusual biotype profiles. If rRNA is high across all
#   samples, repeat library prep with improved rRNA depletion.

ensembl_ids <- rownames(counts_qc)

biotype_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = ensembl_ids,
  keytype = "ENSEMBL",
  columns = c("SYMBOL","GENETYPE")
) %>%
  distinct(ENSEMBL, .keep_all=TRUE) %>%
  mutate(biotype_simple = case_when(
    GENETYPE == "protein-coding" ~ "Protein-coding",
    GENETYPE == "pseudo"         ~ "Pseudogene",
    GENETYPE == "ncRNA"          ~ "ncRNA",
    GENETYPE %in% c("snRNA","snoRNA") ~ "snRNA/snoRNA",
    GENETYPE == "rRNA"           ~ "rRNA",
    GENETYPE == "tRNA"           ~ "tRNA",
    is.na(GENETYPE)              ~ "Unknown",
    TRUE                         ~ "Other"
  ))

counts_biotype <- counts_qc %>%
  as.data.frame() %>%
  rownames_to_column("ENSEMBL") %>%
  left_join(biotype_map %>% dplyr::select(ENSEMBL, biotype_simple),
            by="ENSEMBL") %>%
  mutate(biotype_simple = ifelse(is.na(biotype_simple), "Unknown", biotype_simple)) %>%
  group_by(biotype_simple) %>%
  summarize(across(where(is.numeric), sum)) %>%
  column_to_rownames("biotype_simple")

biotype_pct <- sweep(counts_biotype, 2, colSums(counts_biotype), "/") * 100

biotype_long <- biotype_pct %>%
  rownames_to_column("biotype") %>%
  pivot_longer(-biotype, names_to="sample", values_to="pct") %>%
  left_join(meta_qc %>% dplyr::select(rna_id_fixed, cohort),
            by=c("sample"="rna_id_fixed"))

biotype_order <- biotype_long %>%
  group_by(biotype) %>%
  summarize(mean_pct=mean(pct)) %>%
  arrange(desc(mean_pct)) %>%
  pull(biotype)

biotype_long$biotype <- factor(biotype_long$biotype, levels=rev(biotype_order))

sample_order <- meta_qc %>% arrange(cohort) %>% pull(rna_id_fixed)
biotype_long$sample <- factor(biotype_long$sample, levels=sample_order)

biotype_colors <- c(
  "Protein-coding" = "#2ecc71",
  "ncRNA"          = "#3498db",
  "Pseudogene"     = "#e67e22",
  "snRNA/snoRNA"   = "#9b59b6",
  "rRNA"           = "#e74c3c",
  "tRNA"           = "#f39c12",
  "Other"          = "#95a5a6",
  "Unknown"        = "#bdc3c7"
)

p05 <- ggplot(biotype_long, aes(x=sample, y=pct, fill=biotype)) +
  geom_bar(stat="identity", width=1) +
  geom_hline(yintercept=70, linetype="dashed",
             color="black", linewidth=0.5) +
  annotate("text", x=5, y=72,
           label="70% protein-coding threshold",
           size=3, hjust=0) +
  scale_fill_manual(values=biotype_colors, name="RNA Biotype") +
  labs(
    title    = "05. RNA Biotype Distribution",
    subtitle = paste0("GSE161731 | ", nrow(meta_qc),
                      " samples ordered by cohort | ",
                      "More than 70% protein-coding expected for blood RNA-seq"),
    x = "Sample (ordered by cohort)",
    y = "Percentage of Reads (%)"
  ) +
  theme_minimal(base_size=12) +
  theme(axis.text.x     = element_blank(),
        axis.ticks.x    = element_blank(),
        plot.title      = element_text(face="bold"),
        legend.position = "right")

ggsave("figures/05_biotype_distribution.png", p05,
       width=13, height=6, dpi=300)
cat("Plot 05 saved: biotype distribution\n")


# ==============================================================================
# PLOT 06: PCA — OUTLIER DETECTION
# ==============================================================================
#
# WHAT IT SHOWS:
#   Principal component analysis of all samples using log2(CPM+1) normalized
#   counts, run BEFORE removing outliers. Each point is one sample.
#   Automatic flagging of samples more than 3 SD from the mean on PC1 or PC2.
#
# WHAT GOOD LOOKS LIKE:
#   - Samples cluster loosely by cohort/condition
#   - No single sample far isolated from all others
#   - PC1 + PC2 explain more than 20% of variance
#
# WHAT BAD LOOKS LIKE:
#   - One or more samples far from all others (flagged red)
#     → investigate library size and gene detection for those samples
#   - All samples from one batch cluster away from others
#     → strong batch effect; consider batch correction
#   - Random scatter with no structure → possible labeling errors
#
# ACTION:
#   Samples more than 3 SD from the mean on PC1 or PC2 are candidates for
#   removal. Always cross-reference with library size and gene detection.
#   In this dataset three outliers were identified and removed before
#   downstream analysis (see QC section above).

pct_var <- round(100 * pca_raw$sdev^2 / sum(pca_raw$sdev^2))
pca_df$label <- ifelse(pca_df$outlier, pca_df$sample, "")

p06 <- ggplot(pca_df, aes(x=PC1, y=PC2, color=cohort)) +
  geom_point(size=2.5, alpha=0.8) +
  geom_point(data=pca_df %>% filter(outlier),
             aes(x=PC1, y=PC2), color="red",
             size=4, shape=1, stroke=1.5) +
  geom_text_repel(aes(label=label), size=3,
                  color="red", max.overlaps=20) +
  scale_color_manual(values=palette_cohort, name="Cohort") +
  labs(
    title    = "06. PCA — Outlier Detection",
    subtitle = paste0("GSE161731 | log2(CPM+1) | PC1=", pct_var[1],
                      "%, PC2=", pct_var[2],
                      "% | Red circles = potential outliers (>3 SD from mean)"),
    x = paste0("PC1 (", pct_var[1], "%)"),
    y = paste0("PC2 (", pct_var[2], "%)")
  ) +
  theme_minimal(base_size=12) +
  theme(plot.title = element_text(face="bold"))

ggsave("figures/06_pca_outlier_detection.png", p06,
       width=9, height=6, dpi=300)
cat("Plot 06 saved: PCA outlier detection\n")
cat("Outliers flagged:", sum(pca_df$outlier), "\n")


# ==============================================================================
# SESSION INFO
# ==============================================================================

cat("\n=== Session Info ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("R version:", R.version.string, "\n")
cat("edgeR:", as.character(packageVersion("edgeR")), "\n")
cat("org.Hs.eg.db:", as.character(packageVersion("org.Hs.eg.db")), "\n")
sessionInfo()