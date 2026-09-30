library(tidyverse)

ortho_data <- read_tsv("genes_with_compartments/gill_all_species_TSS_compartments.tsv", show_col_types = FALSE)

ortho_data <- ortho_data |>
  add_count(Orthogroup, spc, name = "n_copies") |>
  filter(Type == "Singleton") |> 
  filter(n_copies == 1)

ortho_data |>
  filter(Type == "Singleton") |>  # adjust column name
  count(Orthogroup, spc, name = "n_copies_per_species") |>
  count(n_copies_per_species)

ortho_data_1to1 <- ortho_data |>
  add_count(Orthogroup, spc, name = "n_copies") |>
  filter(n_copies == 1) |>
  select(-n_copies)

# --- Derive discrete A/B compartment calls from E1 sign ---
# Standard convention: E1 > 0 -> A compartment, E1 < 0 -> B compartment.
# If a `compartment` column already exists upstream, swap this line to use it directly.
ortho_data_1to1 <- ortho_data_1to1 |>
  mutate(compartment = if_else(E1 > 0, "A", "B"))

compartment_wide <- ortho_data_1to1 |>
  select(Orthogroup, spc, compartment) |>
  distinct() |>
  pivot_wider(
    names_from = spc,
    values_from = compartment
  )

compartment_wide |>
  summarise(across(-Orthogroup, ~sum(is.na(.)))) |>
  pivot_longer(everything(), names_to = "species", values_to = "n_missing")


species_list <- setdiff(colnames(compartment_wide), "Orthogroup")

species_pairs <- combn(species_list, 2, simplify = FALSE)

# --- Pairwise % compartment identity (A/B match rate) across singleton orthologs ---
pairwise_conservation <- map_dfr(species_pairs, function(pair) {
  sp1 <- pair[1]
  sp2 <- pair[2]
  
  pair_data <- compartment_wide |>
    select(Orthogroup, all_of(sp1), all_of(sp2)) |>
    drop_na()  # only orthologs with a called compartment in both species
  
  n_pairs <- nrow(pair_data)
  
  if (n_pairs == 0) {
    return(tibble(
      species1 = sp1,
      species2 = sp2,
      n_orthologs = n_pairs,
      n_conserved = NA_integer_,
      pct_conserved = NA_real_,
      p_value = NA_real_
    ))
  }
  
  n_conserved <- sum(pair_data[[sp1]] == pair_data[[sp2]])
  pct_conserved <- 100 * n_conserved / n_pairs
  
  # Test whether conservation exceeds the 50% expected by chance (2 categories)
  binom_test <- binom.test(n_conserved, n_pairs, p = 0.5, alternative = "greater")
  
  tibble(
    species1 = sp1,
    species2 = sp2,
    n_orthologs = n_pairs,
    n_conserved = n_conserved,
    pct_conserved = pct_conserved,
    p_value = binom_test$p.value
  )
})

# --- Multiple testing correction across all pairwise comparisons ---
pairwise_conservation <- pairwise_conservation |>
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    p_adj_display = sprintf("%.10f", p_adj)
  )

# plot
heatmap_data <- pairwise_conservation |>
  select(species1, species2, pct_conserved, n_orthologs) |>
  bind_rows(
    pairwise_conservation |>
      select(species1 = species2, species2 = species1, pct_conserved, n_orthologs)
  )

species_order <- c("Upyg", "Eluc", "Omyk", "Salp", "Ssal")  # adjust to your phylogenetic order

heatmap_data <- heatmap_data |>
  mutate(
    species1 = factor(species1, levels = species_order),
    species2 = factor(species2, levels = species_order)
  )


p1 <- ggplot(heatmap_data, aes(x = species2, y = species1, fill = pct_conserved)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", pct_conserved)),
            color = "black", size = 4.2, fontface = "bold",
            vjust = -0.3) +
  geom_text(aes(label = paste0("n=", n_orthologs)),
            color = "black", size = 2.8, alpha = 0.85,
            vjust = 1.6) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "#DDDDDD", high = "#B40426",
                       midpoint = 50, limits = c(0, 100), name = "% conserved") +
  labs(title = "Singleton ortholog A/B compartment conservation - gill",
       subtitle = "Pairwise % identical compartment call") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(), axis.title = element_blank())
p1
ggsave("2sept_threed/singleton_compartment_conservation_gill.png", p1, width = 8, height = 6, dpi = 200)


library(ape)
species_list <- c("Ssal", "Salp", "Omyk", "Eluc", "Upyg")
outgroup     <- c("Eluc", "Upyg")

# --- 1. Build distance matrix from pairwise % compartment conservation (heatmap_data) ---
# distance = 1 - proportion conserved (0 = identical compartments everywhere, 1 = never matches)

n <- length(species_list)
dist_matrix <- matrix(0, nrow = n, ncol = n,
                      dimnames = list(species_list, species_list))

for (i in seq_len(nrow(heatmap_data))) {
  s1 <- heatmap_data$species1[i]
  s2 <- heatmap_data$species2[i]
  d  <- 1 - (heatmap_data$pct_conserved[i] / 100)
  dist_matrix[s1, s2] <- d
  dist_matrix[s2, s1] <- d
}

compartment_dist <- as.dist(dist_matrix)

# --- 2. Build UPGMA and NJ trees, root NJ on the outgroups ---

upgma_tree <- hclust(compartment_dist, method = "average")
nj_tree    <- nj(compartment_dist)
nj_rooted  <- root(nj_tree, outgroup = outgroup, resolve.root = TRUE)

# --- 3. Prepare the per-ortholog compartment matrix for bootstrapping ---
# compartment_wide: one row per Orthogroup, one column per species (with a leading
# Orthogroup column) -> transpose to species (rows) x orthologs (columns),
# and keep only orthologs with a compartment call across all six species.

compartment_matrix_raw <- as.matrix(compartment_wide[, species_list])
rownames(compartment_matrix_raw) <- compartment_wide$Orthogroup
compartment_matrix <- t(compartment_matrix_raw)
compartment_matrix <- compartment_matrix[, colSums(is.na(compartment_matrix)) == 0]

cat(sprintf("Orthologs used for bootstrap: %d (out of %d total)\n",
            ncol(compartment_matrix), nrow(compartment_wide)))

# --- 4. Bootstrap: resample orthologs, rebuild the tree, repeat ---
# Distance is 1 - pairwise % identical compartment calls, matching step 1.
# Bootstrap runs on the UNROOTED tree (standard, reliable use of
# boot.phylo); each replicate is then rerooted the same way as nj_rooted,
# and support is matched back by CLADE MEMBERSHIP (prop.clades) rather
# than by node position, since rerooting can renumber internal nodes.

build_match_tree <- function(x) {
  # x: species (rows) x orthologs (columns) matrix of "A"/"B" calls
  sp <- rownames(x)
  ns <- length(sp)
  d <- matrix(0, ns, ns, dimnames = list(sp, sp))
  for (i in seq_len(ns - 1)) {
    for (j in (i + 1):ns) {
      pct <- mean(x[i, ] == x[j, ])
      d[i, j] <- 1 - pct
      d[j, i] <- 1 - pct
    }
  }
  nj(as.dist(d))
}

set.seed(1)
boot_result <- boot.phylo(
  nj_tree, compartment_matrix, build_match_tree,
  B = 1000, rooted = FALSE, trees = TRUE
)

boot_trees_rooted <- lapply(
  boot_result$trees,
  function(t) root(t, outgroup = outgroup, resolve.root = TRUE)
)
class(boot_trees_rooted) <- "multiPhylo"

bootstrap_support <- prop.clades(nj_rooted, boot_trees_rooted, rooted = TRUE)
bootstrap_support[is.na(bootstrap_support)] <- 0

# --- 5. Save ONE combined PDF: rooted NJ (with bootstrap) + UPGMA ---

graphics.off()

pdf("2sept_threed/compartment_nj_upgma_trees_gill.pdf", width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(4, 2, 3, 6))

n_tips <- Ntip(nj_rooted)

plot(nj_rooted, main = "Rooted NJ, with bootstrap support",
     cex = 0.9, label.offset = 0.01, no.margin = FALSE,
     y.lim = c(-1, n_tips), font = 4)   # font 4 = bold italic (keeps species-name italics, adds bold)

nodelabels(bootstrap_support, frame = "none", bg = "white",
           adj = c(0.5, -0.8), cex = 0.8, font = 2)   # font 2 = bold

add.scale.bar(x = 0, y = -0.5, cex = 0.8)   # placed in the blank space below all tips

plot(upgma_tree, main = "UPGMA (A/B compartment) - gill",
     xlab = "", sub = "", ylab = "1 - % conserved", cex = 0.9)

par(mfrow = c(1, 1))
dev.off()

cat("Saved: compartment_nj_upgma_trees_gill.pdf\n")

