library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(AnnotationDbi)
library(org.Hs.eg.db)

setwd("~/Research/Data/Fetal_Brain/")

# Setup Output Directories
output <- "~/Research/Figures_and_Results/Fetal_Brain/"
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

#=====================#
# Import RNA-seq Data #
#=====================#
data <- read.delim("Fetal_Brain_gene_counts.txt", header=TRUE, row.names=1, sep="\t", check.names = FALSE)
metadata <- read.delim("fetal_brain_metadata.tsv", header = TRUE, sep = "\t", row.names=1)

# Order data and metadata so sample IDs align
data    <- data[, order(colnames(data))]
metadata <- metadata[order(rownames(metadata)), ]

# Sanity check
stopifnot(all(colnames(data) == rownames(metadata)))

# Filter low coverage Genes
# Compute CPM 
lib_sizes <- colSums(data)
cpm_vals <- t(t(data) / lib_sizes * 1e6)

# Filter: keep genes with CPM > 1 in at least 75% of samples
keep <- rowMeans(cpm_vals > 1) > 0.75
data <- data[keep, ]

 # Barplot summarising number of genes removed with CPM filter
statsbarplot <- tibble(
  Category = factor(c("Before filter", "Removed", "After filter"),
                    levels = c("Before filter", "Removed", "After filter")),
  Genes = c(nrow(cpm_vals), nrow(cpm_vals) - sum(keep), sum(keep))
)

ggplot(statsbarplot, aes(x = Category, y = Genes, fill = Category)) +
  # width = 0.7 makes bars look less blocky; color = "black" adds a crisp outline
  geom_col(width = 0.7, color = "black", linewidth = 0.5) +
  
  # Increased text size, made bold, and slightly shifted up
  geom_text(aes(label = Genes), vjust = -0.8, size = 5, fontface = "bold") +
  
  # Scientific color palette (Okabe-Ito inspired highlight scheme)
  scale_fill_manual(values = c("Before filter" = "#999999",  # Neutral Light Gray
                               "Removed"       = "#D55E00",  # Vermilion (Draws the eye)
                               "After filter"  = "#4D4D4D")) + # Dark Gray
  
  # theme_classic is widely preferred in journals over theme_minimal for bar plots
  theme_classic() + 
  
  # Expand the top of the Y-axis by 15% so the geom_text labels aren't cut off in the PDF
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  
  labs(#title = "Gene Filtering Summary",
       #subtitle = "CPM > 1 in 75% of samples", # Moved the threshold to a subtitle for cleaner look
       x = NULL, 
       y = "Number of genes") +
  
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.text.x = element_text(size = 12, face = "bold", color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)))

ggsave(file.path(output, "Gene_Filtering_Summary.pdf"), width = 8, height = 6, dpi = 1000, device = "pdf")

#==========================#
# PCA Before Normalization #
#==========================#
pca_pre <- prcomp(t(data), center = TRUE, scale. = FALSE)
ve_pre  <- pca_pre$sdev^2 / sum(pca_pre$sdev^2)

pc_pre <- as.data.frame(pca_pre$x) |>
  tibble::rownames_to_column("library_id") |>
  dplyr::left_join(tibble::rownames_to_column(metadata, "library_id"),
                   by = "library_id") |>
  dplyr::mutate(
    Tissue = factor(dplyr::recode(tissue,
                                  "frontal_cortex" = "Frontal cortex",
                                  "spinal_cord"    = "Spinal cord"),
                    levels = c("Frontal cortex", "Spinal cord")),
    Group  = factor(dplyr::recode(group,
                                  "control"            = "Normal morphology",
                                  "neural_tube_defect" = "Neural tube defect"),
                    levels = c("Normal morphology", "Neural tube defect")))

# Okabe–Ito (colour-blind safe)
group_pal <- c("Normal morphology"      = "#0072B2",   # blue
               "Neural tube defect" = "#D55E00")   # vermilion

theme_pub <- theme_classic(base_size = 12) +
  theme(plot.title   = element_text(face = "bold", hjust = 0, size = 13),
        axis.title   = element_text(face = "bold"),
        axis.text    = element_text(colour = "black"),
        legend.title = element_text(face = "bold"),
        legend.position = "right",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        axis.line    = element_blank())

pc_lab <- function(i) sprintf("PC%d  (%.1f%%)", i, 100 * ve_pre[i])

pca_raw <- ggplot(pc_pre, aes(PC1, PC2)) +
  stat_ellipse(aes(group = interaction(Tissue, Group), fill = Group),
               geom = "polygon", alpha = 0.08, colour = NA, level = 0.80,
               show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.3) +
  geom_point(aes(shape = Tissue, fill = Group),
             size = 3.8, stroke = 0.8, colour = "black") +
  ggrepel::geom_text_repel(aes(label = sample_id),
                           size = 2.6, max.overlaps = 12,
                           segment.size = 0.2, colour = "black",
                           show.legend = FALSE) +
  scale_shape_manual(values = c("Frontal cortex" = 21, "Spinal cord" = 24),
                     name   = "Tissue",
                     guide  = guide_legend(
                       order = 1,
                       override.aes = list(fill = "grey60", colour = "black",
                                           size = 3.6))) +
  scale_fill_manual(values = group_pal, name = "Group",
                    guide  = guide_legend(
                      order = 2,
                      override.aes = list(shape = 21, colour = "black",
                                          size = 3.6))) +
  labs(x = pc_lab(1), y = pc_lab(2)) +
  theme_pub

ggsave(file.path(output, "PCA_pre_normalization.pdf"),
       pca_raw, width = 8, height = 6, dpi = 1000, device = cairo_pdf)

#======================#
# DESeq2 Normalization #
#======================#
# 1.  Build the DESeq2 dataset 
# Prepare colData with the experimental factors 
# plot labels will be set to "Normal morphology" / "Neural tube defect"

coldata <- metadata
coldata$Tissue <- factor(dplyr::recode(coldata$tissue,
                                       "frontal_cortex" = "Frontal_cortex",
                                       "spinal_cord"    = "Spinal_cord"),
                         levels = c("Frontal_cortex", "Spinal_cord"))
coldata$Group  <- factor(dplyr::recode(coldata$group,
                                       "control"            = "Normal_morphology",
                                       "neural_tube_defect" = "Neural_tube_defect"),
                         levels = c("Normal_morphology", "Neural_tube_defect"))
coldata$Case   <- factor(coldata$case_id)   # one fetus = one Case; two tissues share it

# Sanity check: counts and metadata must be in the same sample order
stopifnot(all(colnames(data) == rownames(coldata)))

# A single factor encoding Group × Tissue (4 levels)
coldata$Condition <- factor(
  paste(coldata$Group, coldata$Tissue, sep = "_"),
  levels = c("Normal_morphology_Frontal_cortex",
             "Normal_morphology_Spinal_cord",
             "Neural_tube_defect_Frontal_cortex",
             "Neural_tube_defect_Spinal_cord"))

# Build DESeqDataSet, fit model, and normalize 
data <- round(data)
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData   = coldata,
                              design    = ~ Condition)

dds <- DESeq(dds)

# Variance-stabilising transform
# vst() removes the count-mean-variance dependence so the subsequent analysis isn't dominated
# by highly-expressed genes. blind = FALSE lets the transform respect the
# experimental design (recommended for downstream analyses)

vsd <- vst(dds, blind = FALSE)

resultsNames(dds)

#====== NB ==========
# vsd is used only for distance-based downstream applications that are easily skewed by the extreme variance of low-count genes. 
# This includes:
# Principal Component Analysis (PCA) 
# Hierarchical clustering and Heatmaps
# Machine learning models
#=========================# 

#=========================#
# PCA After Normalization #
#=========================#

# PCA on the VST-transformed matrix (all genes, no extra filtering)
pca_post <- prcomp(t(assay(vsd)), center = TRUE, scale. = FALSE)
ve_post  <- pca_post$sdev^2 / sum(pca_post$sdev^2)

pc_post <- as.data.frame(pca_post$x) |>
  tibble::rownames_to_column("library_id") |>
  dplyr::left_join(tibble::rownames_to_column(metadata, "library_id"),
                   by = "library_id") |>
  dplyr::mutate(
    Tissue = factor(dplyr::recode(tissue,
                                  "frontal_cortex" = "Frontal cortex",
                                  "spinal_cord"    = "Spinal cord"),
                    levels = c("Frontal cortex", "Spinal cord")),
    Group  = factor(dplyr::recode(group,
                                  "control"            = "Normal morphology",
                                  "neural_tube_defect" = "Neural tube defect"),
                    levels = c("Normal morphology", "Neural tube defect")))

pc_lab_post <- function(i) sprintf("PC%d  (%.1f%%)", i, 100 * ve_post[i])

pca_norm <- ggplot(pc_post, aes(PC1, PC2)) +
  stat_ellipse(aes(group = interaction(Tissue, Group), fill = Group),
               geom = "polygon", alpha = 0.08, colour = NA, level = 0.80,
               show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.3) +
  geom_point(aes(shape = Tissue, fill = Group),
             size = 3.8, stroke = 0.8, colour = "black") +
  scale_shape_manual(values = c("Frontal cortex" = 21, "Spinal cord" = 24),
                     name   = "Tissue",
                     guide  = guide_legend(
                       order = 1,
                       override.aes = list(fill = "grey60", colour = "black",
                                           size = 3.6))) +
  scale_fill_manual(values = group_pal, name = "Group",
                    guide  = guide_legend(
                      order = 2,
                      override.aes = list(shape = 21, colour = "black",
                                          size = 3.6))) +
  labs(x = pc_lab_post(1), y = pc_lab_post(2)) +
  theme_pub

ggsave(file.path(output, "PCA_post_normalization.pdf"),
       pca_norm, width = 8, height = 6, dpi = 1000, device = cairo_pdf)
#================== End of PCA ================================#

#===================================#
# Remove outlier within spinal cord #
#===================================#
outlier_lib <- "08TF001617_S77"

# Defensive check — fail loudly if the library is already gone
stopifnot(outlier_lib %in% colnames(data))

# Drop from counts, metadata, and coldata
data     <- data[,    colnames(data)    != outlier_lib]
metadata <- metadata[rownames(metadata) != outlier_lib, ]
coldata  <- coldata[ rownames(coldata)  != outlier_lib, ]

# Re-verify alignment
stopifnot(all(colnames(data) == rownames(metadata)))
stopifnot(all(colnames(data) == rownames(coldata)))

# Cohort balance after exclusion
message(sprintf("After exclusion: %d libraries × %d genes",
                ncol(data), nrow(data)))
cat("\nCondition counts:\n");                print(table(coldata$Condition))
cat("\nCase pairing (n libraries per case):\n");  print(table(coldata$Case))

#=========================================#
# Re-build DESeqDataSet on n = 39          #
#=========================================#
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData   = coldata,
                              design    = ~ Condition)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)

resultsNames(dds)

#===========================#
# PCA after outlier removal #
#===========================#
pca_post <- prcomp(t(assay(vsd)), center = TRUE, scale. = FALSE)
ve_post  <- pca_post$sdev^2 / sum(pca_post$sdev^2)

pc_post <- as.data.frame(pca_post$x) |>
  tibble::rownames_to_column("library_id") |>
  dplyr::left_join(tibble::rownames_to_column(metadata, "library_id"),
                   by = "library_id") |>
  dplyr::mutate(
    Tissue = factor(dplyr::recode(tissue,
                                  "frontal_cortex" = "Frontal cortex",
                                  "spinal_cord"    = "Spinal cord"),
                    levels = c("Frontal cortex", "Spinal cord")),
    Group  = factor(dplyr::recode(group,
                                  "control"            = "Normal morphology",
                                  "neural_tube_defect" = "Neural tube defect"),
                    levels = c("Normal morphology", "Neural tube defect")))

pc_lab_post <- function(i) sprintf("PC%d  (%.1f%%)", i, 100 * ve_post[i])

pca_norm <- ggplot(pc_post, aes(PC1, PC2)) +
  stat_ellipse(aes(group = interaction(Tissue, Group), fill = Group),
               geom = "polygon", alpha = 0.08, colour = NA, level = 0.80,
               show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey70", linewidth = 0.3) +
  geom_point(aes(shape = Tissue, fill = Group),
             size = 3.8, stroke = 0.8, colour = "black") +
  scale_shape_manual(values = c("Frontal cortex" = 21, "Spinal cord" = 24),
                     name   = "Tissue",
                     guide  = guide_legend(
                       order = 1,
                       override.aes = list(fill = "grey60", colour = "black",
                                           size = 3.6))) +
  scale_fill_manual(values = group_pal, name = "Group",
                    guide  = guide_legend(
                      order = 2,
                      override.aes = list(shape = 21, colour = "black",
                                          size = 3.6))) +
  labs(x = pc_lab_post(1), y = pc_lab_post(2)) +
  theme_pub

ggsave(file.path(output, "PCA_post_normalization_n39_Outlier_Removed.pdf"),
       pca_norm, width = 8, height = 6, dpi = 1000, device = cairo_pdf)
#==================== End of PCA =======================#

# Data has Ensembl IDs. So we map them onto gene symbols
# Strip the Ensembl version suffix if present (ENSG00000123456.7 -> ENSG00000123456).
# org.Hs.eg.db keys on the unversioned ID
annotate_res <- function(res) {
  df <- as.data.frame(res) |>
    tibble::rownames_to_column("gene_id") |>
    dplyr::arrange(padj)
  
  # Strip the Ensembl version suffix
  df$ensembl_id  <- sub("\\..*$", "", df$gene_id)
  
  # Map using org.Hs.eg.db
  df$gene_symbol <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = df$ensembl_id, column = "SYMBOL",
    keytype = "ENSEMBL", multiVals = "first")
  
  # Manual patch for C4A-AS1, with a fallback to Ensembl ID for any other NAs
  df <- df |>
    dplyr::mutate(
      gene_symbol = dplyr::case_when(
        ensembl_id == "ENSG00000233627" ~ "C4A-AS1",
        is.na(gene_symbol) | gene_symbol == "" ~ ensembl_id,
        TRUE ~ gene_symbol
      )
    )
  
  df
}

# 2. DEG's for the different  contrasts 
# A. NTD vs Normal morphology, WITHIN frontal cortex
res_group_FC <- results(dds, contrast = c("Condition",
  "Neural_tube_defect_Frontal_cortex", "Normal_morphology_Frontal_cortex"),
  alpha = 0.01)

# B. NTD vs Normal morphology, WITHIN spinal cord
res_group_SC <- results(dds, contrast = c("Condition",
  "Neural_tube_defect_Spinal_cord", "Normal_morphology_Spinal_cord"),
  alpha = 0.01)

# C. Spinal cord vs Frontal cortex, WITHIN normal morphology
res_tissue_NM <- results(dds, contrast = c("Condition",
  "Normal_morphology_Spinal_cord", "Normal_morphology_Frontal_cortex"),
  alpha = 0.01)

# D. Spinal cord vs Frontal cortex, WITHIN NTD
res_tissue_NTD <- results(dds, contrast = c("Condition",
  "Neural_tube_defect_Spinal_cord", "Neural_tube_defect_Frontal_cortex"),
  alpha = 0.01)

# E. Overall NTD vs Normal morphology (averaged across both tissues).
# Numeric contrast vector in resultsNames(dds) order — positions are:
#   1 Intercept | 2 Normal_SC | 3 NTD_FC | 4 NTD_SC
res_group_overall <- results(dds, contrast = c(0, -0.5, 0.5, 0.5),
                             alpha = 0.01)

# F. Overall Spinal cord vs Frontal cortex (averaged across both groups)
res_tissue_overall <- results(dds, contrast = c(0,  0.5, -0.5, 0.5),
                              alpha = 0.01)
# Annotate, and summarise DEGs 
contrasts <- list(
  group_FC       = res_group_FC,        # NTD vs Normal in FC
  group_SC       = res_group_SC,        # NTD vs Normal in SC
  tissue_in_NM   = res_tissue_NM,       # SC vs FC in Normal morphology
  tissue_in_NTD  = res_tissue_NTD,      # SC vs FC in NTD
  group_overall  = res_group_overall,   # NTD vs Normal averaged over tissues
  tissue_overall = res_tissue_overall   # SC vs FC averaged over groups
)

contrasts_df <- lapply(contrasts, annotate_res)

# Quick summary: how many genes pass FDR < 0.01 per contrast
n_de <- sapply(contrasts_df, function(d) sum(d$padj < 0.01, na.rm = TRUE))
print(n_de)

# Write each to disk for the supplementary
for (nm in names(contrasts_df)) {
  readr::write_csv(contrasts_df[[nm]],
                   file.path(output, paste0("DE_", nm, ".csv")))
}

# Volcano plot
# Volcano plot for the spinal-cord NTD contrast (group_SC) ---------------
res_df <- contrasts_df$group_SC

volcano_df <- res_df |>
  dplyr::filter(!is.na(padj)) |>
  dplyr::mutate(
    Direction = dplyr::case_when(
      padj < 0.01 & log2FoldChange > 0 ~ "Up in NTD",
      padj < 0.01 & log2FoldChange < 0 ~ "Down in NTD",
      TRUE                              ~ "Not significant"),
    Direction = factor(Direction,
                       levels = c("Up in NTD", "Down in NTD", "Not significant")))

# 27 hits is small — label every one that has a mapped symbol
to_label <- volcano_df |>
  dplyr::filter(Direction != "Not significant",
                !is.na(gene_symbol), gene_symbol != "") |>
  dplyr::arrange(padj)

direction_pal <- c("Up in NTD"       = "#D55E00",
                   "Down in NTD"     = "#0072B2",
                   "Not significant" = "grey78")

p_volcano <- ggplot(volcano_df, aes(log2FoldChange, -log10(padj))) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey80", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_point(aes(colour = Direction,
                 size   = Direction != "Not significant",
                 alpha  = Direction != "Not significant")) +
  ggrepel::geom_text_repel(data = to_label,
                           aes(label = gene_symbol),
                           size = 2.9, max.overlaps = Inf,
                           segment.size = 0.25, segment.colour = "grey50",
                           min.segment.length = 0,
                           box.padding = 0.5, point.padding = 0.2,
                           colour = "black", show.legend = FALSE) +
  scale_colour_manual(values = direction_pal, name = NULL) +
  scale_size_manual( values = c(`FALSE` = 1.2, `TRUE` = 2.2), guide = "none") +
  scale_alpha_manual(values = c(`FALSE` = 0.45, `TRUE` = 0.95), guide = "none") +
  labs(#title = "NTD vs Normal morphology — spinal cord",
       x = expression(log[2]~fold~change~"("*NTD~vs.~Normal~morphology*")"),
       y = expression(-log[10]~italic(P)[adj])) +
  theme_pub

ggsave(file.path(output, "Volcano_NTD_vs_NM_spinal_cord.pdf"),
       p_volcano, width = 7.5, height = 6, dpi = 1000, device = cairo_pdf)
