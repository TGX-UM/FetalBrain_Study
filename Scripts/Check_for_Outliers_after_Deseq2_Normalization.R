#=========================================#
# Sample-to-sample Pearson correlation    #
#=========================================#
# Within-tissue mean correlation is the principal outlier metric in a
# two-tissue cohort.  A library that correlates with its same-tissue
# siblings at < 0.85 is correlating with them about as poorly as it
# correlates with the opposite tissue.


V       <- assay(vsd)
cor_mat <- cor(V, method = "pearson")

# Build a samp table joined to tissue (uses `metadata`, not `coldata`)
samp <- tibble::tibble(
  library_id = colnames(V),
  sample_id  = metadata[colnames(V), "sample_id"],
  tissue     = metadata[colnames(V), "tissue"]
) |>
  dplyr::mutate(
    Tissue = factor(dplyr::recode(tissue,
                                  "frontal_cortex" = "Frontal cortex",
                                  "spinal_cord"    = "Spinal cord"),
                    levels = c("Frontal cortex", "Spinal cord")))

# Per-library mean correlation, split into within-tissue and across-tissue
mean_within  <- vapply(samp$library_id, function(id) {
  tt <- samp$tissue[samp$library_id == id]
  partners <- samp$library_id[samp$tissue == tt & samp$library_id != id]
  mean(cor_mat[id, partners])
}, numeric(1))
mean_across <- vapply(samp$library_id, function(id) {
  tt <- samp$tissue[samp$library_id == id]
  partners <- samp$library_id[samp$tissue != tt]
  mean(cor_mat[id, partners])
}, numeric(1))

samp <- samp |>
  dplyr::mutate(
    mean_r_within  = mean_within,
    mean_r_across  = mean_across,
    # robust z-score (within tissue) for cohort-adaptive flagging
    z_within = ave(mean_r_within, Tissue,
                   FUN = function(x) (x - median(x)) / mad(x)),
    flag_threshold = mean_r_within < 0.85,
    flag_robust_z  = z_within < -3,
    suspect        = flag_threshold | flag_robust_z)

# Cohort summary per tissue
samp |>
  dplyr::group_by(Tissue) |>
  dplyr::summarise(median_within_r = median(mean_r_within),
                   min_within_r    = min(mean_r_within),
                   median_across_r = median(mean_r_across),
                   n_below_0.85    = sum(flag_threshold),
                   n_robust_z_fail = sum(flag_robust_z)) |>
  print()

# Suspect list, sorted by within-tissue r
suspect_tbl <- samp |>
  dplyr::filter(suspect) |>
  dplyr::arrange(mean_r_within) |>
  dplyr::select(library_id, sample_id, Tissue,
                mean_r_within, mean_r_across, z_within,
                flag_threshold, flag_robust_z)
message(sprintf("\nLibraries failing at least one criterion: %d / %d",
                nrow(suspect_tbl), nrow(samp)))
print(suspect_tbl)


# Bar plot & Heatmap of Outlier

#--- bar chart: within-tissue mean r per library, with 0.85 line -------------
samp_ord <- samp |>
  dplyr::arrange(Tissue, mean_r_within) |>
  dplyr::mutate(library_id = factor(library_id, levels = library_id))

tissue_pal <- c("Frontal cortex" = "#56B4E9",   # sky blue
                "Spinal cord"    = "#E69F00")   # orange

p_bar <- ggplot(samp_ord,
                aes(x = library_id, y = mean_r_within, fill = Tissue)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0.85, linetype = "dashed",
             colour = "#D55E00", linewidth = 0.5) +
  annotate("text", x = nrow(samp_ord), y = 0.86,
           label = "r = 0.85", hjust = 1, vjust = 0,
           colour = "#D55E00", fontface = "bold", size = 3.4) +
  ggrepel::geom_text_repel(data = dplyr::filter(samp_ord, suspect),
                           aes(label = sample_id), size = 2.8,
                           min.segment.length = 0, segment.size = 0.2,
                           nudge_y = 0.015, colour = "black") +
  scale_fill_manual(values = tissue_pal, name = "Tissue") +
  coord_cartesian(ylim = c(min(0.70, min(samp_ord$mean_r_within) - 0.02), 1)) +
  labs(x = NULL, y = "Mean Pearson r vs same-tissue libraries") +
  theme_pub +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())

ggsave(file.path(output, "Sample_correlation_within_tissue.pdf"),
       p_bar, width = 9, height = 6, dpi = 1000, device = cairo_pdf)

#--- annotated heatmap: order rows/cols by tissue → group → library_id -------
order_ids <- samp |>
  dplyr::left_join(tibble::rownames_to_column(metadata, "library_id") |>
                     dplyr::select(library_id, group),
                   by = "library_id") |>
  dplyr::arrange(Tissue, group, library_id) |>
  dplyr::pull(library_id)

n_fc          <- sum(samp$Tissue == "Frontal cortex")
divider_x     <- n_fc + 0.5
divider_y     <- nrow(samp) - n_fc + 0.5

heat_df <- as.data.frame(cor_mat[order_ids, order_ids]) |>
  tibble::rownames_to_column("y_lib") |>
  tidyr::pivot_longer(-y_lib, names_to = "x_lib", values_to = "r") |>
  dplyr::mutate(x_lib = factor(x_lib, levels = order_ids),
                y_lib = factor(y_lib, levels = rev(order_ids)))

p_heat <- ggplot(heat_df, aes(x_lib, y_lib, fill = r)) +
  geom_tile() +
  geom_vline(xintercept = divider_x, colour = "white", linewidth = 0.6) +
  geom_hline(yintercept = divider_y, colour = "white", linewidth = 0.6) +
  scale_fill_viridis_c(option = "viridis", name = "Pearson r",
                       limits = c(min(heat_df$r), 1)) +
  coord_equal() +
  labs(x = NULL, y = NULL) +
  theme_pub +
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
        axis.text.y  = element_text(size = 5),
        panel.border = element_blank(),
        axis.line    = element_blank())

ggsave(file.path(output, "Sample_correlation_heatmap.pdf"),
       p_heat, width = 10, height = 9, dpi = 1000, device = cairo_pdf)
