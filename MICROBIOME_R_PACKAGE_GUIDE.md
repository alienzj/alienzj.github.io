# Comprehensive R Package Guide for Microbiome & Bioinformatics Analysis

**Version:** 1.0
**Date:** 2026-03-20
**Target Audience:** Computational biologists, microbiome researchers, bioinformaticians

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Core Data Structures](#2-core-data-structures)
3. [Data Import & Processing](#3-data-import--processing)
4. [Diversity Analysis](#4-diversity-analysis)
5. [Ordination & Visualization](#5-ordination--visualization)
6. [Differential Abundance Testing](#6-differential-abundance-testing)
7. [Multivariable Association](#7-multivariable-association)
8. [Machine Learning & Classification](#8-machine-learning--classification)
9. [Network Analysis](#9-network-analysis)
10. [Compositional Data Analysis](#10-compositional-data-analysis)
11. [Simulation & Power Analysis](#11-simulation--power-analysis)
12. [Future Directions](#12-future-directions)

---

## 1. Introduction

### 1.1 The Microbiome R Ecosystem

The R ecosystem for microbiome analysis has evolved into a comprehensive framework addressing the unique challenges of microbiome data:

**Key Challenges:**
- **Compositionality**: Relative abundance data (sums to 1 or 100%)
- **Sparsity**: Many zero counts (zero-inflation)
- **High-dimensional**: Thousands of taxa, few samples
- **Phylogenetic structure**: Evolutionary relationships matter
- **Hierarchical taxonomy**: Kingdom → Phylum → Class → Order → Family → Genus → Species
- **Multi-omics integration**: Metagenomics, metatranscriptomics, metabolomics

**Ecosystem Overview:**

```
┌─────────────────────────────────────────────────────────────────┐
│                    Microbiome R Ecosystem                        │
├─────────────────────────────────────────────────────────────────┤
│  Bioconductor Core                                              │
│  ├─ phyloseq (data structure)                                  │
│  ├─ DESeq2 (differential abundance)                            │
│  ├─ edgeR (differential abundance)                             │
│  ├─ Maaslin2 (multivariable association)                       │
│  ├─ SIAMCAT (ML classification)                                │
│  └─ ALDEx2 (compositional analysis)                            │
├─────────────────────────────────────────────────────────────────┤
│  CRAN Packages                                                  │
│  ├─ vegan (community ecology)                                  │
│  ├─ microbiome (tidy workflow)                                 │
│  └─ MicrobiotaProcess (tidy framework)                         │
├─────────────────────────────────────────────────────────────────┤
│  Specialized Tools                                              │
│  ├─ picante (phylogenetic ecology)                             │
│  ├─ iNEXT (diversity estimation)                               │
│  ├─ SpiecEasi (network inference)                              │
│  └─ compositions (compositional data)                          │
└─────────────────────────────────────────────────────────────────┘
```

### 1.2 Choosing the Right Package

| Analysis Goal | Recommended Packages |
|--------------|---------------------|
| Data import & storage | phyloseq, MicrobiotaProcess |
| Alpha diversity | vegan, iNEXT, microbiome |
| Beta diversity | vegan, phyloseq |
| Differential abundance | Maaslin2, DESeq2, ALDEx2, ANCOM-BC, edgeR |
| Multi-variable analysis | Maaslin2 |
| Classification/ML | SIAMCAT |
| Network analysis | SpiecEasi, microbiome |
| Compositional analysis | ALDEx2, zCompositions, compositions |
| Phylogenetic analysis | picante, phyloseq |

---

## 2. Core Data Structures

### 2.1 phyloseq: The Standard Data Container

**Problem Solved:**
Microbiome data comes from multiple sources (OTU tables, taxonomy, metadata, phylogeny) that need to be kept synchronized. phyloseq provides a unified S4 class that ensures data integrity.

**Design Philosophy:**
- **S4 Object System**: Enforces strict class definitions
- **Component Integration**: Links OTU table, taxonomy, metadata, phylogeny
- **Bioconductor Compatibility**: Works with SummarizedExperiment ecosystem
- **Method Dispatch**: Generic functions work across components

**Data Structure:**

```r
library(phyloseq)

# phyloseq object contains 4-5 components:
# 1. otu_table: Matrix of counts (features × samples)
# 2. tax_table: Matrix of taxonomic classification
# 3. sample_data: Data frame of sample metadata
# 4. phylo_tree: Phylogenetic tree (optional)
# 5. refseq: Sequence data (optional)

ps <- phyloseq(otu_table_obj, tax_table_obj, sample_data_obj, phylo_tree_obj)
```

**Key Functions:**

```r
# Access components
otu_table(ps)        # Get OTU table
tax_table(ps)        # Get taxonomy
sample_data(ps)      # Get metadata
phylo_tree(ps)       # Get phylogenetic tree
refseq(ps)           # Get sequences

# Modify components
otu_table(ps) <- new_otu
sample_data(ps) <- new_meta

# Subset
ps_filtered <- subset_samples(ps, Disease == "Case")
ps_taxa <- subset_taxa(ps, Kingdom == "Bacteria")

# Merge
ps_merged <- merge_phyloseq(ps1, ps2)
```

### 2.2 MPSE (MicrobiotaProcess): Tidy Framework

**Problem Solved:**
Traditional microbiome workflows are fragmented. MPSE provides a tidyverse-compatible interface that integrates with modern R workflows.

**Design Philosophy:**
- **Tidyverse Integration**: Pipe-friendly, dplyr-compatible
- **MPSE Class**: MicrobiomeProcessSample class for unified data handling
- **Workflow Abstraction**: Encapsulates common analysis steps

**Usage:**

```r
library(MicrobiotaProcess)

# Create MPSE object
mpse_obj <- MPSE(
    count_data = count_matrix,
    tax_data = tax_matrix,
    sample_data = metadata
)

# Extract distance matrix
dist_matrix <- mpse_obj %>% mp_extract_dist(distmethod = "bray")

# Transform data
mpse_transformed <- mpse_obj %>% mp_transform(method = "clr")
```

### 2.3 SIAMCAT Object

**Problem Solved:**
Machine learning with microbiome data requires careful handling of abundance data, metadata, and labels. SIAMCAT provides a structured object for ML workflows.

**Design Philosophy:**
- **Label-Centric**: Organizes data around prediction targets
- **Feature Selection**: Built-in feature importance methods
- **Cross-Validation**: Integrated validation workflows
- **Confounder Control**: Handles covariates in ML

**Structure:**

```r
siamcat_obj <- list(
    abund_data = abundance_matrix,
    metadata = metadata,
    label = "Disease",
    features = selected_features,
    training_results = list(),
    metrics = list()
)
```

---

## 3. Data Import & Processing

### 3.1 Importing Data

**phyloseq Import Functions:**

```r
# From BIOM format (common from QIIME, mothur)
ps <- import_biom("otu_table.biom")

# From QIIME output
ps <- import_qiime(
    otu_table = "otu_table.txt",
    taxonomy = "taxonomy.txt",
    sample_data = "sample_data.txt",
    suppress_warnings = TRUE
)

# From mothur
ps <- import_mothur(
    shared_file = "shared_file.txt",
    taxonomy_file = "taxonomy.txt"
)

# From general matrices
ps <- phyloseq(
    otu_table(otu_matrix, taxa_are_rows = TRUE),
    tax_table(tax_matrix),
    sample_data(meta_data)
)
```

**Data Format Requirements:**

```r
# OTU Table: Features × Samples or Samples × Features
# Must have unique row/column names
# Values must be non-negative integers

# Taxonomy Table: Features × Taxonomic Ranks
# Rows must match OTU table features
# Format: Kingdom;Phylum;Class;Order;Family;Genus;Species

# Sample Data: Samples × Covariates
# Row names must match OTU table columns
# Can include continuous and categorical variables
```

### 3.2 Data Processing & Transformation

**Common Transformations:**

```r
library(phyloseq)

# 1. Filter low-abundance taxa
ps_filtered <- filter_taxa(
    ps,
    function(x) sum(x > 0) > 10,  # Keep taxa present in >10 samples
    TRUE
)

# 2. Agglomerate at higher taxonomic levels
ps_glom <- tax_glom(ps_filtered, taxrank = "Phylum")

# 3. Transform counts
ps_log <- transform_sample_counts(ps, function(x) log10(x + 1))
ps_sqrt <- transform_sample_counts(ps, sqrt)
ps_prop <- transform_sample_counts(ps, function(x) x / sum(x))

# 4. Rarefaction (subsampling)
ps_rarefied <- rarefy_even_depth(ps, sample.size = 10000)

# 5. CLR transformation (compositional)
ps_clr <- transform_sample_counts(ps, function(x) {
    clr(x + 1e-6)  # Add pseudocount before CLR
})
```

**MicrobiotaProcess Transformations:**

```r
# Multiple transformation methods
mp_transform(mpse_obj, method = "clr")    # Centered log-ratio
mp_transform(mpse_obj, method = "log")     # Log transformation
mp_transform(mpse_obj, method = "sqrt")    # Square root
mp_transform(mpse_obj, method = "hellinger") # Hellinger
```

---

## 4. Diversity Analysis

### 4.1 Alpha Diversity

**Problem Solved:**
Alpha diversity measures the diversity within a single sample. Key questions: How many species? How evenly distributed?

**vegan Implementation:**

```r
library(vegan)

# Extract community matrix
comm_matrix <- otu_table(ps)
if (taxa_are_rows(ps)) {
    comm_matrix <- t(comm_matrix)  # Samples × Features
}

# 1. Species Richness (S)
richness <- specnumber(comm_matrix)

# 2. Shannon Diversity
shannon <- diversity(comm_matrix, index = "shannon")

# 3. Simpson Diversity
simpson <- diversity(comm_matrix, index = "simpson")

# 4. Inverse Simpson
inv_simpson <- diversity(comm_matrix, index = "invsimpson")

# 5. Berger-Parker Index
berger_parker <- diversity(comm_matrix, index = "bergerparker")

# 6. Pielou Evenness
pielou <- diversity(comm_matrix, index = "pielou")

# Multiple metrics at once
richness_metrics <- estimate_richness(
    comm_matrix,
    measures = c("Shannon", "Simpson", "Chao1", "ACE", "Jack1")
)
```

**phyloseq Integration:**

```r
# Direct from phyloseq object
alpha_div <- estimate_richness(ps, measures = c("Shannon", "Simpson", "Chao1"))

# Add to sample data
sample_data(ps)$Shannon <- alpha_div$Shannon
sample_data(ps)$Simpson <- alpha_div$Simpson

# Statistical comparison
library(ggplot2)
ggplot(sample_data(ps), aes(x = Disease, y = Shannon)) +
    geom_boxplot() +
    geom_jitter(width = 0.2) +
    stat_compare_means()  # From ggpubr
```

**iNEXT: Interpolation/Extrapolation**

```r
library(iNEXT)

# Prepare data
sample_list <- lapply(split(otu_table(ps), col(otu_table(ps))), as.vector)

# iNEXT analysis
out <- iNEXT(
    sample_list,
    datatype = "abundance",
    size = c(100:1000, 1000:5000),  # Sample sizes
    nboot = 200                     # Bootstrap replicates
)

# Plot diversity curves
plot(out, type = 1)  # Sample-based
plot(out, type = 2)  # Distance-based

# Extract diversity estimates
summary(out)
```

### 4.2 Beta Diversity

**Problem Solved:**
Beta diversity measures the difference in community composition between samples. Key questions: Are groups different? What drives differences?

**Distance Metrics:**

```r
# Bray-Curtis (abundance-based)
dist_bray <- vegdist(comm_matrix, method = "bray")

# Jaccard (presence/absence)
dist_jaccard <- vegdist(comm_matrix, method = "jaccard")

# UniFrac (phylogenetic)
dist_unifrac_weighted <- UniFrac(ps, weighted = TRUE)
dist_unifrac_unweighted <- UniFrac(ps, weighted = FALSE)

# Canberra distance
dist_canberra <- vegdist(comm_matrix, method = "canberra")

# Euclidean distance
dist_euclidean <- vegdist(comm_matrix, method = "euclidean")

# Hellinger distance
comm_hell <- decostand(comm_matrix, method = "hellinger")
dist_hell <- vegdist(comm_hell, method = "euclidean")
```

**phyloseq Distance Functions:**

```r
# Calculate distances
dist_bray <- distance(ps, method = "bray")
dist_jaccard <- distance(ps, method = "jaccard")
dist_unifrac <- UniFrac(ps, weighted = FALSE, normalized = TRUE)

# Access distance matrix
as.matrix(dist_bray)
```

### 4.3 Statistical Testing of Beta Diversity

**PERMANOVA (adonis2):**

```r
library(vegan)

# Basic PERMANOVA
adonis_result <- adonis2(
    comm_matrix ~ Disease + Age + BMI,
    data = metadata,
    method = "bray",
    permutations = 999,
    by = "margin"  # Marginal effects
)

# With phylogenetic distance
adonis_unifrac <- adonis2(
    dist_unifrac ~ Disease,
    data = metadata,
    permutations = 999
)

# Output format
#               Df SumOfSqs      R2    F.Model        Pr(>F)
# Disease        1   0.2456 0.12345   12.345      0.001 **
# Age            1   0.1234 0.06789    6.123      0.023 *
# Residuals    197   3.9876 0.80866
```

**ANOSIM (Analysis of Similarities):**

```r
# ANOSIM
anosim_result <- anosim(
    dist_bray,
    metadata$Disease,
    permutations = 999
)

# R statistic: 0 = no difference, 1 = complete separation
anosim_result$statistic
anosim_result$signif
```

**PERMDISP (Beta Diversity Dispersion):**

```r
# Test for homogeneity of multivariate dispersions
beta_disp <- betadisper(dist_bray, metadata$Disease)

# Test significance
permutest(beta_disp)

# Plot
plot(beta_disp)
boxplot(beta_disp)

# Important: PERMANOVA can be significant due to dispersion differences
# Always check PERMDISP to interpret PERMANOVA correctly
```

---

## 5. Ordination & Visualization

### 5.1 Ordination Methods

**NMDS (Non-metric Multidimensional Scaling):**

```r
library(vegan)

# Basic NMDS
nmds <- metaMDS(
    comm_matrix,
    distance = "bray",
    k = 2,              # Number of dimensions
    trymax = 20,        # Maximum trials
    autotransform = FALSE,
    wascores = TRUE
)

# Extract coordinates
nmds_coords <- scores(nmds)

# Stress value (lower is better, <0.2 acceptable)
nmds$stress

# Plot
plot(nmds, type = "n")
points(nmds, display = "sites", col = metadata$Disease)
text(nmds, display = "species", col = "gray")
```

**PCoA (Principal Coordinates Analysis):**

```r
# From phyloseq
ord <- ordinate(ps, method = "PCoA", distance = "bray")

# From vegan
eigenvals <- cmdscale(dist_bray, k = 10, eig = TRUE)
pcoa_coords <- eigenvals$points
```

**RDA (Redundancy Analysis):**

```r
# Constrained ordination
rda_result <- rda(
    comm_matrix ~ Disease + Age + pH,
    data = metadata,
    scale = TRUE
)

# Variance partitioning
varpart_result <- varpart(
    comm_matrix,
    ~ Disease + Age,
    ~ pH + Temperature
)
plot(varpart_result)

# Forward selection
rda_forward <- ordiR2step(
    rda(comm_matrix ~ ., data = metadata),
    scope = formula(rda(comm_matrix ~ 1, data = metadata)),
    direction = "forward",
    R2scope = TRUE
)
```

### 5.2 Visualization with ggplot2

**Taxonomic Bar Plot:**

```r
library(phyloseq)
library(ggplot2)

# Using plot_bar
p <- plot_bar(ps, fill = "Phylum") +
    facet_wrap(~Disease) +
    theme_bw() +
    theme(legend.position = "right")

# Customized with ggplot2
ps_melt <- psmelt(ps)

ggplot(ps_melt, aes(x = Sample, y = Abundance, fill = Phylum)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_y_continuous(labels = scales::percent) +
    facet_wrap(~Disease) +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```

**Ordination Plot:**

```r
# Using phyloseq
plot_ordination(ps, ord, color = "Disease", shape = "Gender") +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(aes(group = Disease), type = "t", linetype = 2) +
    theme_minimal()

# Using vegan + ggplot2
nmds_df <- as.data.frame(nmds$points)
nmds_df$Disease <- metadata$Disease

ggplot(nmds_df, aes(x = MDS1, y = MDS2, color = Disease)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_text(aes(label = Sample), hjust = 0, vjust = 0) +
    theme_minimal()

# Add environmental vectors
env_fit <- envfit(nmds ~ pH + Temperature, data = metadata, permutations = 999)
plot(env_fit, p.max = 0.05, col = "red")
```

**Heatmap:**

```r
# Using phyloseq
plot_heatmap(ps, x = "Disease", sort = "groups")

# Custom heatmap with pheatmap
library(pheatmap)

# Subset top taxa
top_taxa <- names(sort(rowSums(otu_table(ps)), decreasing = TRUE))[1:50]
heatmap_data <- t(log10(otu_table(ps)[top_taxa, ] + 1))

pheatmap(
    heatmap_data,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = metadata,
    scale = "row"
)
```

---

## 6. Differential Abundance Testing

### 6.1 DESeq2

**Problem Solved:**
Differential abundance testing with count data following negative binomial distribution. Handles library size differences and biological variability.

**Design Philosophy:**
- **Negative Binomial GLM**: Models count data appropriately
- **Shrinkage Estimation**: Stabilizes dispersion and fold change estimates
- **Multiple Testing Correction**: Built-in FDR control
- **Flexible Design**: Supports complex experimental designs

**Full Implementation:**

```r
library(DESeq2)

# 1. Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
    countData = otu_table(ps),      # Integer count matrix
    colData = sample_data(ps),      # Metadata
    design = ~ Disease              # Formula
)

# 2. Pre-filtering (optional but recommended)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# 3. Run DESeq pipeline
dds <- DESeq(
    dds,
    fitType = "parametric",         # "parametric", "local", "mean", "glmGamPoi"
    betaPrior = FALSE,
    parallel = FALSE,
    minReplicatesForReplace = 7
)

# 4. Extract results
res <- results(
    dds,
    contrast = c("Disease", "Case", "Control"),
    lfcThreshold = 0,
    altHypothesis = "greaterAbs",
    pAdjustMethod = "BH",
    alpha = 0.05
)

# 5. Shrink log2 fold changes
res_shrunk <- lfcShrink(
    dds,
    coef = "Disease_Case_vs_Control",
    type = "apeglm",               # "apeglm", "ashr", "normal"
    lfcThreshold = 0
)

# 6. Results table
results_table <- as.data.frame(res_shrunk)
sig_features <- results_table[results_table$padj < 0.05 & abs(results_table$log2FoldChange) > 1, ]

# 7. Visualization
# MA plot
plotMA(res, main = "DESeq2 Results", alpha = 0.05)

# PCA
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "Disease")

# Heatmap of top features
top_features <- head(order(res$padj), 50)
pheatmap(
    assay(vsd)[top_features, ],
    annotation_col = sample_data(dds)
)
```

**Complex Designs:**

```r
# Multi-factor design
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ Batch + Disease + Age
)

# Interaction term
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ Treatment * Time
)

# Likelihood ratio test for multiple factors
dds_full <- DESeqDataSetFromMatrix(counts, metadata, ~ Treatment + Time + Treatment:Time)
dds_reduced <- DESeqDataSetFromMatrix(counts, metadata, ~ Treatment + Time)
dds <- DESeq(dds_full, test = "LRT", reduced = ~ Treatment + Time)
```

### 6.2 edgeR

**Problem Solved:**
Differential expression/abundance analysis for count data. Alternative to DESeq2 with different dispersion estimation methods.

**Design Philosophy:**
- **TMM Normalization**: Trimmed Mean of M-values
- **Empirical Bayes**: Shrinkage of dispersion estimates
- **Flexible Testing**: Exact test (2 groups) or GLM (multiple groups)

**Full Implementation:**

```r
library(edgeR)

# 1. Create DGEList
y <- DGEList(
    counts = otu_table(ps),
    group = metadata$Disease
)

# 2. Filter lowly expressed features
y <- filterByExpr(y)

# 3. Calculate normalization factors
y <- calcNormFactors(
    y,
    method = "TMM",               # "TMM", "TMMW", "RLE", "upperquartile"
    logratioTrim = 0.3,
    sumTrim = 0.05
)

# 4. Design matrix
design <- model.matrix(~Disease, data = metadata)

# 5. Estimate dispersion
y <- estimateDisp(y, design, robust = TRUE)

# 6. Exact test (two groups)
et <- exactTest(
    y,
    pair = c("Control", "Case"),
    prior.count = 0.125
)

# 7. Results
top_tags <- topTags(et, n = Inf, adjust.method = "BH")
sig_features <- top_tags$table[top_tags$table$FDR < 0.05, ]

# GLM approach (multiple factors)
fit <- glmQLFit(y, design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = 2)
top_tags_glm <- topTags(qlf, n = Inf)

# Visualization
# MD plot
plotMD(et, column = 1)

# Heatmap
top_features <- order(et$table$PValue)[1:50]
cpm_values <- cpm(y)[top_features, ]
pheatmap(cpm_values)
```

### 6.3 ALDEx2

**Problem Solved:**
Compositional data analysis using centered log-ratio transformation with Monte Carlo sampling. Designed specifically for microbiome data.

**Design Philosophy:**
- **Compositional Awareness**: CLR transformation accounts for relative nature
- **Monte Carlo**: Multiple imputation of technical variation
- **Robust Statistics**: Welch's t-test or Wilcoxon rank-sum
- **Effect Sizes**: Reports magnitude of differences

**Full Implementation:**

```r
library(ALDEx2)

# 1. Basic differential abundance
x <- aldex(
    reads = otu_table(ps),        # Integer count matrix (features × samples)
    conditions = metadata$Disease, # Character vector of group labels
    mc.samples = 128,             # Monte Carlo iterations (more = more accurate)
    denom = "all",                # "all" or specific taxon name
    test = "t",                   # "t", "wilcox", "glm", "kr"
    eff = TRUE,                   # Calculate effect sizes
    alt.hypothesis = "two.sided", # "two.sided", "greater", "less"
    verbose = FALSE
)

# 2. Get results
aldex_results <- aldex.ttest(x)

# Results columns:
# - effect: Effect size (mean difference in CLR space)
# - we.ep: Welch's p-value
# - wi.eBH: Wilcoxon BH-corrected p-value
# - median, mean: Abundance in each group

sig_features <- aldex_results[aldex_results$we.eBH < 0.05 & abs(aldex_results$effect) > 0.5, ]

# 3. With design matrix (adjust for covariates)
design <- model.matrix(~Disease + Age + BMI, data = metadata)
x_glm <- aldex(
    reads = otu_table(ps),
    design = design,
    mc.samples = 128
)

glm_results <- aldex.glm(x_glm, formula = ~Disease)

# 4. Correlation analysis
corr_results <- aldex.corr(x, method = "spearman", mc.reps = 100)

# 5. Visualization
# Effect size plot
effect_size <- aldex_results[aldex_results$we.eBH < 0.05, ]
ggplot(effect_size, aes(x = reorder(feature, effect), y = effect)) +
    geom_point() +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed")
```

### 6.4 ANCOM-BC

**Problem Solved:**
Differential abundance with bias correction for compositional data. Controls false discovery rate while accounting for sampling fraction differences.

**Design Philosophy:**
- **Bias Correction**: Estimates and corrects for sampling fraction bias
- **Zero Handling**: Models structural zeros
- **Compositional**: Works in log-ratio space
- **Robust FDR Control**: Holm or FDR adjustment

**Full Implementation:**

```r
library(ANCOMBC)

# 1. Using phyloseq object
out <- ancombc(
    phyloseq = ps,
    formula = "Disease",              # Main effect
    zero_cut = 0.9,                   # Max proportion of zeros to keep
    lib_cut = 0,                      # Min library size
    group = NULL,                     # Stratification variable
    corr_method = "holm",             # P-value correction
    th_logfc = 0,                     # LogFC threshold
    adj_pmethod = "holm",             # Adj p-value method
    pv_adj_method = "fdr",            # P-value adjustment
    dmm = FALSE,                      # Dirichlet-multinomial model
    random_effect = NULL,             # Random effects
    stratified = NULL,                # Stratified analysis
    logfc_threshold = 0,              # Minimum logFC
    verbose = TRUE
)

# 2. Results
res <- out$res
# Columns:
# - diff_abn: Differential abundance (TRUE/FALSE)
# - W: Test statistic
# - p_val: Raw p-value
# - q_val: Adjusted p-value
# - log2FC: Log2 fold change
# - avg: Mean abundance

sig_taxa <- res[res$q_val < 0.05 & res$diff_abn == TRUE, ]

# 3. Normalized feature table
normalized_table <- out$feature_table

# 4. Visualization
# Volcano plot
ggplot(res, aes(x = log2FC, y = -log10(q_val))) +
    geom_point(alpha = 0.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed")
```

### 6.5 Comparison of Differential Abundance Methods

| Method | Best For | Pros | Cons |
|--------|----------|------|------|
| **DESeq2** | Count data, robust | Well-established, handles complex designs | May be conservative |
| **edgeR** | Small sample sizes | Flexible, good power | Similar to DESeq2 |
| **ALDEx2** | Compositional data | CLR-based, robust to zeros | Computationally intensive |
| **ANCOM-BC** | FDR control | Bias correction, good specificity | May be conservative |
| **Maaslin2** | Multi-variable | Handles covariates, random effects | Requires larger samples |

---

## 7. Multivariable Association

### 7.1 MaAsLin2

**Problem Solved:**
Identify microbial features associated with metadata while controlling for confounders. The go-to tool for clinical microbiome studies.

**Design Philosophy:**
- **General Linear Models**: LM, MLM, GLM, GAM, SVM, RF
- **Normalization**: Multiple options (TSS, CLR, etc.)
- **Random Effects**: Handle repeated measures
- **Multiple Testing**: FDR control built-in
- **Transformations**: LOG, AITCHISON, NONE

**Full Implementation:**

```r
library(Maaslin2)

# 1. Prepare data
# Count matrix: Features × Samples
# Metadata: Samples × Covariates

# 2. Run MaAsLin2
results <- Maaslin2(
    input_data = otu_table(ps),       # Abundance matrix
    input_metadata = sample_data(ps), # Metadata
    output = "maaslin2_output",       # Output directory
    formula = "Disease + Age + BMI",  # Fixed effects
    random_effect = "PatientID",      # Random effects (optional)
    normalization = "TSS",            # "TSS", "TSSp", "CLR", "TSSp"
    transform = "LOG",                # "NONE", "LOG", "AITCHISON", "LOGIT"
    analysis_method = "LM",           # "LM", "MLM", "LR", "SVM", "RF", "GBM"
    q_value_correction = "BH",        # P-value adjustment
    min_abundance = 0,                # Min abundance threshold
    min_prevalence = 0,               # Min prevalence threshold
    max_significance = 0.1,           # Max q-value for significance
    zero_cut = 0.9,                   # Max proportion of zeros
    out_overdispersion = TRUE,        # Output overdispersion
    verbose = TRUE
)

# 3. Results structure
# results$results: Main results table
# results$random: Random effect results
# results$metadata: Input metadata

# 4. Results table columns:
# - feature: Feature name
# - variable: Covariate
# - coef: Coefficient (effect size)
# - std_error: Standard error
# - t_stat: T-statistic
# - p_val: Raw p-value
# - q_val: Adjusted p-value
# - intercept: Intercept

# 5. Filter significant associations
sig_results <- results$results[results$results$q_val < 0.05, ]

# 6. Extract specific associations
disease_associations <- sig_results[sig_results$variable == "Disease", ]

# 7. Visualization
# Association plot
plot_maaslin_results(
    results,
    q_value_threshold = 0.1,
    min_abs_coef = 0.5
)
```

**Advanced Usage:**

```r
# Multiple testing with different methods
results_fdr <- Maaslin2(
    input_data = otu_table(ps),
    input_metadata = sample_data(ps),
    formula = "Disease",
    q_value_correction = "fdr"
)

results_bh <- Maaslin2(
    input_data = otu_table(ps),
    input_metadata = sample_data(ps),
    formula = "Disease",
    q_value_correction = "BH"
)

# Random forest analysis
results_rf <- Maaslin2(
    input_data = otu_table(ps),
    input_metadata = sample_data(ps),
    formula = "Disease",
    analysis_method = "RF",
    n_estimators = 1000
)
```

---

## 8. Machine Learning & Classification

### 8.1 SIAMCAT

**Problem Solved:**
Machine learning for microbiome-based prediction and feature selection. Handles high-dimensional sparse data with proper validation.

**Design Philosophy:**
- **Label-Centric**: Organizes around prediction target
- **Multiple ML Methods**: GLM, glmnet, SVM, Random Forest
- **Feature Selection**: Built-in methods for selecting predictive taxa
- **Confounder Adjustment**: Control for covariates
- **Comprehensive Metrics**: AUC, accuracy, F1, etc.

**Full Implementation:**

```r
library(SIAMCAT)

# 1. Create SIAMCAT object
siamcat_obj <- create_siamcat(
    abund_data = otu_table(ps),         # Features × Samples
    metadata = sample_data(ps),         # Sample metadata
    label = "Disease"                   # Column name for prediction
)

# 2. Preprocess data
siamcat_obj <- add_metadata(
    siamcat_obj,
    metadata = sample_data(ps),
    label = "Disease"
)

# 3. Feature filtering
siamcat_obj <- filter_features(
    siamcat_obj,
    min_prevalence = 0.1,              # Min prevalence
    min_abundance = 1e-6               # Min abundance
)

# 4. Train model
siamcat_obj <- train_siamcat(
    siamcat_obj,
    method = "glmnet",                  # "glm", "glmnet", "svm", "rf"
    label = "Disease",
    test_id = "Disease_Classifier",
    cores = 4,
    cv_folds = 10,                      # Cross-validation folds
    alpha = 1,                          # LASSO (1) or Ridge (0)
    nlambda = 100                       # Number of lambda values
)

# 5. Calculate metrics
siamcat_obj <- calc_metrics(
    siamcat_obj,
    metric = "auc",                     # "auc", "accuracy", "f1", "aucpr"
    test_id = "Disease_Classifier"
)

# 6. Feature selection
siamcat_obj <- select_features(
    siamcat_obj,
    method = "LR",                      # "LR", "RF", "MDS"
    n_features = 20,
    test_id = "Disease_Classifier"
)

# 7. Holdout testing
siamcat_obj <- holdout_test(
    siamcat_obj,
    label = "Disease",
    test_prop = 0.2,
    n_repeats = 10
)

# 8. Visualization
# ROC curve
plotROC(siamcat_obj, test_id = "Disease_Classifier")

# Feature importance
plotFeatureImportance(
    siamcat_obj,
    test_id = "Disease_Classifier",
    n_features = 20
)

# Confusion matrix
plotConfusionMatrix(siamcat_obj, test_id = "Disease_Classifier")
```

**Multi-class Classification:**

```r
# For multi-class labels
siamcat_obj <- create_siamcat(
    abund_data = otu_table(ps),
    metadata = sample_data(ps),
    label = "Disease_Subtype"  # Multiple disease subtypes
)

# One-vs-rest classification
siamcat_obj <- train_siamcat(
    siamcat_obj,
    method = "glmnet",
    label = "Disease_Subtype",
    one_vs_rest = TRUE
)
```

---

## 9. Network Analysis

### 9.1 SpiecEasi

**Problem Solved:**
Infer microbial association networks from composition data. Distinguishes direct from indirect associations.

**Design Philosophy:**
- **Sparse Inverse Covariance**: Graphical lasso for network inference
- **Compositional**: Handles relative abundance data
- **Multiple Methods**: Meinshausen-Bühlmann, neighborhood selection
- **Stability Selection**: Robust edge inference

**Implementation:**

```r
library(SpiecEasi)

# 1. Prepare data
# CLR transformation recommended
comm_matrix <- otu_table(ps)
comm_clr <- t(apply(comm_matrix, 1, function(x) {
    log(x / exp(mean(log(x + 1e-10)))) + 1e-10
}))

# 2. Fit network
spiec_easi_obj <- spiec.easi(
    comm_clr,
    method = "mb",                  # "mb", "glasso", "nealasso"
    lambda.max.ratio = 0.01,
    nlambda = 20,
    ncirc = 100,
    pulsar.p = 0.05,                # P-value threshold
    usePulsar = TRUE,
    check = TRUE,
    parallel = TRUE
)

# 3. Extract network
network <- getRefit(spiec_easi_obj)

# 4. Stability selection
pulsar_results <- getPulsar(spiec_easi_obj)

# 5. Visualization
plot(spiec_easi_obj)

# Export network
write.table(network, "network.txt", sep = "\t")
```

---

## 10. Compositional Data Analysis

### 10.1 zCompositions

**Problem Solved:**
Zero replacement for compositional data analysis. Microbiome data has many zeros that prevent log-transformation.

**Implementation:**

```r
library(zCompositions)

# 1. Multiple imputation for zeros
data_zc <- cmultRepl(
    data = otu_table(ps),
    method = "CZM",                 # "CZM", "IM", "LRM", "QRN", "SM"
    label = "zeros",
    seed = 1234
)

# 2. Single imputation
data_single <- cmultRepl(
    otu_table(ps),
    method = "QRN",
    label = "zeros"
)

# 3. Perform CLR on imputed data
clr_data <- clr(data_zc$data)
```

### 10.2 compositions

**Problem Solved:**
General compositional data analysis tools for R.

**Implementation:**

```r
library(compositions)

# 1. CLR transformation
clr_data <- clr(comm_matrix + 1e-6)

# 2. Aitchison distance
ait_dist <- aitchison(comm_matrix)

# 3. Co-variation analysis
covar_data <- coVar(clr_data)
```

---

## 11. Simulation & Power Analysis

### 11.1 SIMCAT (Simulation)

**Problem Solved:**
Generate synthetic microbiome datasets for method testing and power analysis.

**Note:** SIMCAT may refer to different packages. The main simulation tools are:
- **microbiomeSim**: Simulate microbiome data
- **micSim**: Microbiome simulation

**General Simulation Approach:**

```r
# Simulate with realistic properties
# - Zero inflation
# - Overdispersion
# - Taxonomic correlations
# - Phylogenetic structure

# Use packages like microbiomeSim or generate synthetic data manually
```

---

## 12. Future Directions

### 12.1 Current Trends

**1. Tidyverse Integration:**
```r
# Move toward pipe-friendly APIs
ps %>%
    filter_taxa(...) %>%
    transform(...) %>%
    estimate_richness()
```

**2. Single-Cell Microbiome:**
- Integration with SingleCellExperiment
- scMetagenomics methods

**3. Multi-omics Integration:**
- mixOmics integration
- Multi-table analysis

**4. Deep Learning:**
- Microbiome-specific neural networks
- Autoencoders for dimensionality reduction

**5. Cloud-Native Analysis:**
- Integration with bioconductor cloud
- Parallel computing improvements

### 12.2 Emerging Packages to Watch

| Package | Status | Focus Area |
|---------|--------|------------|
| **microbiomeSeq** | Active | Tidy workflow |
| **microme** | Development | Modern microbiome analysis |
| **gmic** | Emerging | Graph-based microbiome analysis |
| **PhyloFactor** | Active | Phylogenetic partitioning |
| **corncob** | Active | Beta-binomial regression |

### 12.3 Best Practices Evolution

```r
# Future direction: More compositional-aware methods
# Current: CLR + standard tests
# Future: Direct compositional models (robCompositions, etc.)

# Better integration with longitudinal data
# Current: Maaslin2 random effects
# Future: Dedicated longitudinal microbiome models

# Explainable AI for microbiome ML
# Current: Feature importance
# Future: SHAP values, LIME for microbiome
```

---

## Appendix: Quick Reference

### Installation

```r
# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "phyloseq",
    "DESeq2",
    "edgeR",
    "Maaslin2",
    "SIAMCAT",
    "ALDEx2",
    "ANCOMBC",
    "MicrobiotaProcess",
    "SpiecEasi",
    "zCompositions"
))

# CRAN packages
install.packages(c(
    "vegan",
    "ggplot2",
    "dplyr",
    "tidyr",
    "pheatmap"
))
```

### Complete Analysis Workflow

```r
# Load libraries
library(phyloseq)
library(vegan)
library(Maaslin2)
library(ggplot2)

# 1. Import data
ps <- import_biom("data.biom")

# 2. Filter and transform
ps <- filter_taxa(ps, rowSums(otu_table(ps)) > 10, TRUE)
ps_log <- transform_sample_counts(ps, function(x) log10(x + 1))

# 3. Diversity analysis
alpha_div <- estimate_richness(ps)
dist_bray <- distance(ps, "bray")

# 4. Statistical testing
adonis2(dist_bray ~ Disease, data = sample_data(ps))

# 5. Differential abundance
results <- Maaslin2(
    input_data = otu_table(ps),
    input_metadata = sample_data(ps),
    formula = "Disease",
    output = "output"
)

# 6. Visualization
plot_bar(ps, fill = "Phylum")
plot_ordination(ps, ordinate(ps, "PCoA", "bray"), color = "Disease")
```

---

*This documentation is continuously updated. For the latest versions and features, check the respective package vignettes and documentation.*
