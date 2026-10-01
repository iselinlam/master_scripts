library(tidyverse)
library(ape)

species_list <- c("Ssal", "Salp", "Omyk", "Eluc", "Upyg")
outgroup     <- c("Eluc", "Upyg")

spc_levels    <- c("Tthy", "Salp", "Ssal", "Omyk")

species_names <- tribble(
  ~spc, ~common_name, ~scientific, 
  "Ssal", "Atlantic Salmon", "Salmo salar",
  "Salp", "Artic Char", "Salvelinus_alpinus",
  "Omyk", "Rainbow Trout", "Oncorhynchus mykiss",
  "Tthy", "European Grayling", "Thymallus thymallus", 
  "Eluc", "Northern Pike", "Esox lucius",
  "Upyg", "Eastern Mudminnow", "Umbra pygmaea"
)

output_dir <- "singleton_analysis"
dir.create(output_dir, showWarnings = FALSE)

# comp_dir      <- "genes_assigned_compartment_zscore"
tissue_levels <- c("liver", "brain", "gill")

# Reading input files and filter to only 1to1 singletons

# make_E1_z_wide <- function(tissue) {
#   
#   files <- list.files(comp_dir,
#                       pattern = paste0("_genes_wTSS_compartment_", tissue, "\\.tsv$"),
#                       full.names = TRUE)
#   
#   names(files) <- sub("_genes_wTSS_compartment_.*$", "", basename(files))
#   
#   message(tissue, ": ", length(files), " species files (",
#           paste(names(files), collapse = ", "), ")")
#   
#   ortho_data <- map_dfr(files,
#                         ~ read_tsv(.x, col_types = cols(.default = col_character())),
#                         .id = "spc") |>
#     mutate(E1_z = as.numeric(E1_z))
#   
#   singletons <- ortho_data |>
#     filter(Type == "Singleton") |>
#     add_count(Orthogroup, spc, name = "n_copies")
#   
#   # count table: all singleton genes vs 1:1 singleton genes, per species
#   counts <- singletons |>
#     group_by(spc) |>
#     summarise(n_singleton      = n(),
#               n_1to1_singleton = sum(n_copies == 1),
#               pct_1to1         = round(100 * n_1to1_singleton / n_singleton, 1),
#               .groups = "drop") |>
#     mutate(tissue = tissue, .before = 1)
#   
#   ortho_data_1to1 <- singletons |>
#     filter(n_copies == 1) |>
#     select(-n_copies)
#   
#   E1_z_wide <- ortho_data_1to1 |>
#     select(Orthogroup, spc, E1_z) |>
#     distinct() |>
#     pivot_wider(names_from = spc, values_from = E1_z)
#   
#   cat("\n--", tissue, "missingness --\n")
#   E1_z_wide |>
#     summarise(across(-Orthogroup, ~ sum(is.na(.)))) |>
#     pivot_longer(everything(), names_to = "species", values_to = "n_missing") |>
#     print(n = Inf)
#   
#   write_tsv(E1_z_wide,
#             file.path(output_dir, paste0("E1_z_singletons_1to1_", tissue, ".tsv")))
#   
#   list(E1_z_wide = E1_z_wide, counts = counts)
# }
# 
# res <- set_names(tissue_levels) |> map(make_E1_z_wide)
# 
# E1_z_wide_list <- map(res, "E1_z_wide")   # E1_z_wide_list$gill etc.
# 
# singleton_counts <- map_dfr(res, "counts") |>
#   mutate(tissue = factor(tissue, levels = tissue_levels)) |>
#   arrange(tissue, spc)
# 
# print(singleton_counts, n = Inf)
# write_tsv(singleton_counts, file.path(output_dir, "E1_z_singleton_1to1_counts.tsv"))


# ---- input files: 
# - E1_z_singletons_1to1_brain.tsv
# - E1_z_singletons_1to1_gill.tsv
# - E1_z_singletons_1to1_liver.tsv

# Pairwise correlation and Spearman heatmap, per tissue

species_order_all <- c("Upyg", "Eluc", "Tthy", "Omyk", "Salp", "Ssal")

plot_E1_z_cor <- function(tissue) {
  
  E1_z_wide <- read_tsv(file.path(output_dir, paste0("E1_z_singletons_1to1_", tissue, ".tsv")),
                      show_col_types = FALSE)
  
  species_order <- intersect(species_order_all, setdiff(names(E1_z_wide), "Orthogroup"))
  species_pairs <- combn(species_order, 2, simplify = FALSE)
  
  pairwise_cor <- map_dfr(species_pairs, function(pair) {
    sp1 <- pair[1]
    sp2 <- pair[2]
    
    pair_data <- E1_z_wide |>
      select(Orthogroup, all_of(sp1), all_of(sp2)) |>
      drop_na()
    
    n_pairs <- nrow(pair_data)
    if (n_pairs < 3) {
      return(tibble(species1 = sp1, species2 = sp2, n_orthologs = n_pairs,
                    pearson_r = NA_real_, spearman_rho = NA_real_,
                    p_value = NA_real_, p_value_spearman = NA_real_))
    }
    
    pearson_test  <- cor.test(pair_data[[sp1]], pair_data[[sp2]], method = "pearson")
    spearman_test <- cor.test(pair_data[[sp1]], pair_data[[sp2]],
                              method = "spearman", exact = FALSE)
    
    tibble(
      species1 = sp1, species2 = sp2, n_orthologs = n_pairs,
      pearson_r = unname(pearson_test$estimate),
      spearman_rho = unname(spearman_test$estimate),
      p_value = pearson_test$p.value,
      p_value_spearman = spearman_test$p.value
    )
  }) |>
    mutate(p_adj = p.adjust(p_value, method = "BH"),
           p_adj_spearman = p.adjust(p_value_spearman, method = "BH"),
           tissue = tissue, .before = 1)
  
  write_tsv(pairwise_cor,
            file.path(output_dir, paste0("E1_z_pairwise_cor_singleton_1to1_", tissue, ".tsv")))
  
  # mirror to a full symmetric table for the heatmap
  heatmap_data <- bind_rows(
    pairwise_cor |> select(species1, species2, pearson_r, spearman_rho, n_orthologs),
    pairwise_cor |> select(species1 = species2, species2 = species1,
                           pearson_r, spearman_rho, n_orthologs)
  ) |>
    mutate(species1 = factor(species1, levels = species_order),
           species2 = factor(species2, levels = species_order))
  
  p <- ggplot(heatmap_data, aes(x = species2, y = species1, fill = spearman_rho)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.2f", spearman_rho)),
              color = "black", size = 4.2, fontface = "bold", vjust = -0.3) +
    geom_text(aes(label = paste0("n=", n_orthologs)),
              color = "black", size = 2.8, alpha = 0.85, vjust = 1.6) +
    scale_fill_gradient2(low = "#3B4CC0", mid = "#DDDDDD", high = "#B40426",
                         midpoint = 0.5, limits = c(0, 1),
                         oob = scales::squish, name = "spearman_rho") +
    labs(title = paste("Singleton ortholog E1_z conservation -", tissue),
         subtitle = "Pairwise Spearman correlation") +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(), axis.title = element_blank())
  
  print(p)
  
  ggsave(file.path(output_dir, paste0("E1_z_singleton_conservation_", tissue, ".png")),
         p, width = 8, height = 6, dpi = 200)
  
  list(cor = pairwise_cor, plot = p)
}

cor_res <- set_names(tissue_levels) |> map(plot_E1_z_cor)

# all tissues in one table
all_pairwise_cor <- map_dfr(cor_res, "cor")
write_tsv(all_pairwise_cor, file.path(output_dir, "E1_z_pairwise_cor_singleton_1to1_all_tissues.tsv"))


# Distance matrix (1 - Spearman rho) and NJ + UPGMA and bootstrap, per tissue

# tip label lookup: "Atlantic Salmon (Ssal)"
tip_labels <- setNames(paste0(species_names$common_name, " (", species_names$spc, ")"),
                       species_names$spc)

make_trees <- function(tissue, B = 1000) {
  
  pw <- read_tsv(file.path(output_dir, paste0("E1_z_pairwise_cor_singleton_1to1_", tissue, ".tsv")),
                 show_col_types = FALSE)
  
  sp_present <- intersect(species_order_all, unique(c(pw$species1, pw$species2)))
  outgrp     <- intersect(outgroup, sp_present)
  
  # --- distance matrix (1 - Spearman rho) ---
  dist_matrix <- matrix(0, nrow = length(sp_present), ncol = length(sp_present),
                        dimnames = list(sp_present, sp_present))
  for (i in seq_len(nrow(pw))) {
    d <- 1 - pw$spearman_rho[i]
    dist_matrix[pw$species1[i], pw$species2[i]] <- d
    dist_matrix[pw$species2[i], pw$species1[i]] <- d
  }
  write_tsv(as_tibble(dist_matrix, rownames = "species"),
            file.path(output_dir, paste0("E1_z_dist_matrix_", tissue, ".tsv")))
  
  E1_z_dist    <- as.dist(dist_matrix)
  upgma_tree <- hclust(E1_z_dist, method = "average")
  nj_tree    <- nj(E1_z_dist)
  
  # --- unrooted diagnostic ---
  nj_unrooted <- nj_tree
  nj_unrooted$tip.label <- unname(tip_labels[nj_unrooted$tip.label])
  
  graphics.off()
  pdf(file.path(output_dir, paste0("E1_z_nj_unrooted_", tissue, ".pdf")), width = 7, height = 6)
  plot(nj_unrooted, type = "unrooted",
       main = paste("Unrooted NJ (check Eluc/Upyg placement) -", tissue))
  dev.off()
  
  # --- root on outgroup; fall back to Upyg alone if not monophyletic ---
  nj_rooted <- tryCatch(
    root(nj_tree, outgroup = outgrp, resolve.root = TRUE),
    error = function(e) {
      cat("NOTE (", tissue, "): outgroup not monophyletic, rooting on Upyg alone.\n")
      root(nj_tree, outgroup = "Upyg", resolve.root = TRUE)
    }
  )
  
  # --- bootstrap on orthologs ---
  E1_z_wide <- read_tsv(file.path(output_dir, paste0("E1_z_singletons_1to1_", tissue, ".tsv")),
                      show_col_types = FALSE)
  E1_z_matrix <- t(as.matrix(E1_z_wide[, sp_present]))
  colnames(E1_z_matrix) <- E1_z_wide$Orthogroup
  
  cat(sprintf("%s: orthologs available for bootstrap: %d\n", tissue, ncol(E1_z_matrix)))
  
  build_corr_tree <- function(x) {
    cor_matrix <- cor(t(x), use = "pairwise.complete.obs", method = "spearman")
    nj(as.dist(1 - cor_matrix))
  }
  
  set.seed(1)
  boot_result <- boot.phylo(nj_tree, E1_z_matrix, build_corr_tree,
                            B = B, rooted = FALSE, trees = TRUE)
  
  boot_trees_rooted <- lapply(
    boot_result$trees,
    function(t) tryCatch(root(t, outgroup = outgrp, resolve.root = TRUE),
                         error = function(e) NULL)
  )
  n_failed <- sum(vapply(boot_trees_rooted, is.null, logical(1)))
  if (n_failed > 0) {
    cat(sprintf("%s: %d of %d bootstrap replicates excluded (outgroup not monophyletic)\n",
                tissue, n_failed, length(boot_trees_rooted)))
  }
  boot_trees_rooted <- boot_trees_rooted[!vapply(boot_trees_rooted, is.null, logical(1))]
  class(boot_trees_rooted) <- "multiPhylo"
  
  bootstrap_support <- prop.clades(nj_rooted, boot_trees_rooted, rooted = TRUE)
  bootstrap_support[is.na(bootstrap_support)] <- 0
  
  # --- relabel tips with common names (after bootstrap matching) ---
  nj_rooted$tip.label <- unname(tip_labels[nj_rooted$tip.label])
  upgma_tree$labels   <- unname(tip_labels[upgma_tree$labels])
  
  # --- save rooted NJ (bootstrap) + UPGMA ---
  graphics.off()
  pdf(file.path(output_dir, paste0("E1_z_upgma_", tissue, ".pdf")), width = 12, height = 7)
  par(mfrow = c(1, 2), mar = c(4, 2, 3, 6))
  
  plot(nj_rooted, main = paste("Rooted NJ, with bootstrap support -", tissue),
       cex = 0.9, label.offset = 0.01, no.margin = FALSE,
       y.lim = c(-1, Ntip(nj_rooted)), font = 4)
  nodelabels(bootstrap_support, frame = "none", adj = c(0.5, -0.8),
             cex = 0.8, font = 2)
  add.scale.bar(x = 0, y = -0.5, cex = 0.8)
  
  plot(upgma_tree, main = paste("UPGMA -", tissue),
       xlab = "", sub = "", ylab = "1 - Spearman rho", cex = 0.9)
  
  par(mfrow = c(1, 1))
  dev.off()
  
  list(dist = dist_matrix, nj = nj_rooted, upgma = upgma_tree,
       support = bootstrap_support)
}

tree_res <- set_names(tissue_levels) |> map(make_trees)

cat("Saved to", output_dir, "/ for each tissue:\n")
cat("  E1_z_nj_unrooted_<tissue>.pdf\n")
cat("  E1_z_upgma_<tissue>.pdf  (rooted NJ with bootstrap + UPGMA)\n")
cat("  E1_z_dist_matrix_<tissue>.tsv\n")
