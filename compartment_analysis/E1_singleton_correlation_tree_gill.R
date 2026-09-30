# =============================================================================
# Singleton ortholog E1 compartment conservation — gill
# Data prep -> pairwise correlation heatmap -> UPGMA/NJ distance trees ->
# bootstrapped rooted NJ. All outputs saved to singleton_analysis/, with a
# "_gill" suffix on every file.
# =============================================================================

library(tidyverse)
library(ape)

output_dir <- "singleton_analysis"
dir.create(output_dir, showWarnings = FALSE)

species_list <- c("Ssal", "Salp", "Omyk", "Eluc", "Upyg")
outgroup     <- c("Eluc", "Upyg")

# =============================================================================
# --- 1. Load data, restrict to 1:1 singleton orthologs ---
# =============================================================================

ortho_data <- read_tsv("genes_with_compartments/gill_all_species_TSS_compartments.tsv",
                       show_col_types = FALSE)

ortho_data_1to1 <- ortho_data |>
  filter(Type == "Singleton") |>
  add_count(Orthogroup, spc, name = "n_copies") |>
  filter(n_copies == 1) |>
  select(-n_copies)

# --- E1 matrix: one row per Orthogroup, one column per species ---

e1_wide <- ortho_data_1to1 |>
  select(Orthogroup, spc, E1) |>
  distinct() |>
  pivot_wider(names_from = spc, values_from = E1)

# quick missingness check per species
e1_wide |>
  summarise(across(-Orthogroup, ~ sum(is.na(.)))) |>
  pivot_longer(everything(), names_to = "species", values_to = "n_missing") |>
  print(n = Inf)

# =============================================================================
# --- 2. Pairwise Pearson + Spearman correlation across all species pairs ---
# =============================================================================

species_pairs <- combn(species_list, 2, simplify = FALSE)

pairwise_cor <- map_dfr(species_pairs, function(pair) {
  sp1 <- pair[1]
  sp2 <- pair[2]
  
  pair_data <- e1_wide |>
    select(Orthogroup, all_of(sp1), all_of(sp2)) |>
    drop_na()
  
  n_pairs <- nrow(pair_data)
  if (n_pairs < 3) {
    return(tibble(species1 = sp1, species2 = sp2, n_orthologs = n_pairs,
                  pearson_r = NA_real_, spearman_rho = NA_real_, p_value = NA_real_))
  }
  
  pearson_test  <- cor.test(pair_data[[sp1]], pair_data[[sp2]], method = "pearson")
  spearman_test <- cor.test(pair_data[[sp1]], pair_data[[sp2]], method = "spearman")
  
  tibble(
    species1 = sp1, species2 = sp2, n_orthologs = n_pairs,
    pearson_r = pearson_test$estimate,
    spearman_rho = spearman_test$estimate,
    p_value = pearson_test$p.value
  )
})

pairwise_cor <- pairwise_cor |>
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# mirror to a full symmetric table for the heatmap
heatmap_data <- bind_rows(
  pairwise_cor |> select(species1, species2, pearson_r, spearman_rho, n_orthologs),
  pairwise_cor |> select(species1 = species2, species2 = species1,
                         pearson_r, spearman_rho, n_orthologs)
)

species_order <- c("Upyg", "Eluc", "Omyk", "Salp", "Ssal")
heatmap_data <- heatmap_data |>
  mutate(species1 = factor(species1, levels = species_order),
         species2 = factor(species2, levels = species_order))

# =============================================================================
# --- 3. Heatmap ---
# =============================================================================

p1 <- ggplot(heatmap_data, aes(x = species2, y = species1, fill = spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", spearman_rho)),
            color = "black", size = 4.2, fontface = "bold", vjust = -0.3) +
  geom_text(aes(label = paste0("n=", n_orthologs)),
            color = "black", size = 2.8, alpha = 0.85, vjust = 1.6) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "#DDDDDD", high = "#B40426",
                       midpoint = 0.5, limits = c(0, 1), name = "spearman_rho") +
  labs(title = "Singleton ortholog E1 conservation - gill",
       subtitle = "Pairwise Spearman correlation") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(), axis.title = element_blank())
p1

ggsave(file.path(output_dir, "singleton_E1_conservation_gill.png"),
       p1, width = 8, height = 6, dpi = 200)

# =============================================================================
# --- 4. Distance matrix (1 - Spearman rho) and trees ---
# =============================================================================

n <- length(species_list)
dist_matrix <- matrix(0, nrow = n, ncol = n,
                      dimnames = list(species_list, species_list))

for (i in seq_len(nrow(heatmap_data))) {
  s1 <- heatmap_data$species1[i]
  s2 <- heatmap_data$species2[i]
  d  <- 1 - heatmap_data$spearman_rho[i]
  dist_matrix[s1, s2] <- d
  dist_matrix[s2, s1] <- d
}

e1_dist <- as.dist(dist_matrix)

upgma_tree <- hclust(e1_dist, method = "average")
nj_tree    <- nj(e1_dist)

# --- unrooted diagnostic (check Eluc/Upyg group together before rooting) ---

graphics.off()
pdf(file.path(output_dir, "e1_spearman_nj_unrooted_gill.pdf"), width = 7, height = 6)
plot(nj_tree, type = "unrooted", main = "Unrooted NJ (check Eluc/Upyg placement) - gill")
dev.off()

# --- root on Eluc+Upyg; fall back to Eluc alone if they aren't monophyletic ---

nj_rooted <- tryCatch(
  root(nj_tree, outgroup = outgroup, resolve.root = TRUE),
  error = function(e) {
    cat("NOTE: Eluc+Upyg are not monophyletic in this tree — falling back to",
        "rooting on Upyg alone (Eluc has been seen nesting WITH an ingroup",
        "species, e.g. Tthy, in some runs — check the unrooted plot for this",
        "specific dataset before trusting this fallback).\n")
    root(nj_tree, outgroup = "Upyg", resolve.root = TRUE)
  }
)

# =============================================================================
# --- 5. Bootstrap: resample orthologs, rebuild the tree, repeat ---
# Uses pairwise-complete Spearman correlation per replicate, matching how
# heatmap_data / the point-estimate tree were built (no complete-case-only
# restriction, which would needlessly shrink the resampling pool).
# =============================================================================

e1_matrix_raw <- as.matrix(e1_wide[, species_list])
rownames(e1_matrix_raw) <- e1_wide$Orthogroup
e1_matrix <- t(e1_matrix_raw)

cat(sprintf("Orthologs available for bootstrap resampling: %d\n", ncol(e1_matrix)))

build_corr_tree <- function(x) {
  cor_matrix <- cor(t(x), use = "pairwise.complete.obs", method = "spearman")
  nj(as.dist(1 - cor_matrix))
}

set.seed(1)
boot_result <- boot.phylo(
  nj_tree, e1_matrix, build_corr_tree,
  B = 1000, rooted = FALSE, trees = TRUE
)

# reroot each replicate the same way; drop any that aren't monophyletic
boot_trees_rooted <- lapply(
  boot_result$trees,
  function(t) tryCatch(root(t, outgroup = outgroup, resolve.root = TRUE),
                       error = function(e) NULL)
)

n_failed <- sum(vapply(boot_trees_rooted, is.null, logical(1)))
if (n_failed > 0) {
  cat(sprintf("%d of %d bootstrap replicates excluded (Eluc+Upyg not monophyletic)\n",
              n_failed, length(boot_trees_rooted)))
}
boot_trees_rooted <- boot_trees_rooted[!vapply(boot_trees_rooted, is.null, logical(1))]
class(boot_trees_rooted) <- "multiPhylo"

# match support to nj_rooted's nodes by clade membership, robust to rerooting
bootstrap_support <- prop.clades(nj_rooted, boot_trees_rooted, rooted = TRUE)
bootstrap_support[is.na(bootstrap_support)] <- 0

# =============================================================================
# --- 6. Save rooted NJ (with bootstrap) + UPGMA together ---
# =============================================================================

graphics.off()
pdf(file.path(output_dir, "e1_spearman_nj_rooted_upgma_gill.pdf"), width = 12, height = 7)
par(mfrow = c(1, 2), mar = c(4, 2, 3, 6))

n_tips <- Ntip(nj_rooted)

plot(nj_rooted, main = "Rooted NJ, with bootstrap support - gill",
     cex = 0.9, label.offset = 0.01, no.margin = FALSE,
     y.lim = c(-1, n_tips), font = 4)
nodelabels(bootstrap_support, frame = "circle", bg = "white",
           adj = c(0.5, -0.8), cex = 0.8, font = 2)
add.scale.bar(x = 0, y = -0.5, cex = 0.8)

plot(upgma_tree, main = "UPGMA - gill",
     xlab = "", sub = "", ylab = "1 - Spearman rho", cex = 0.9)

par(mfrow = c(1, 1))
dev.off()

cat("Saved to singleton_analysis/:\n")
cat("  singleton_E1_conservation_gill.png\n")
cat("  e1_spearman_nj_unrooted_gill.pdf\n")
cat("  e1_spearman_nj_rooted_upgma_gill.pdf\n")

