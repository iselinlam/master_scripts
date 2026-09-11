library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(readr)


df <- read_tsv("genes_with_compartments/brain_all_species_TSS_compartments.tsv", show_col_types = FALSE)

species_list     <- c("Ssal", "Salp", "Omyk", "Tthy", "Eluc", "Upyg")
ss4r_species     <- c("Ssal", "Salp", "Omyk", "Tthy")
outgroup_species <- c("Eluc", "Upyg")

# --- 1. Restrict to clean type-per-species-class rows (as before) ---

raw <- df[
  (df$spc %in% ss4r_species     & df$Type == "SS4Ronly") |
    (df$spc %in% outgroup_species & df$Type == "Singleton"),
  c("HOG", "spc", "compartment")
]

# --- 2. Collapse to one consensus call per HOG x species ---

consensus_fun <- function(x) {
  u <- unique(x)
  if (length(u) == 1) u else "Divergent"
}

consensus_long <- aggregate(compartment ~ HOG + spc, data = raw, FUN = consensus_fun)

# --- 3. Reshape to one row per HOG, one column per species ---

wide <- reshape(consensus_long, idvar = "HOG", timevar = "spc", direction = "wide")
names(wide) <- gsub("^compartment\\.", "", names(wide))

# --- 4. Diagnostics: how often do the two ohnolog copies disagree? ---
# (only meaningful for SS4R species — outgroups have 1 copy, never "Divergent")

cat("Ohnolog disagreement rate per SS4R species (Divergent / total HOGs with data):\n")
for (s in ss4r_species) {
  n_total     <- sum(!is.na(wide[[s]]))
  n_divergent <- sum(wide[[s]] == "Divergent", na.rm = TRUE)
  cat(sprintf("  %s: %d / %d (%.1f%%)\n", s, n_divergent, n_total, 100 * n_divergent / n_total))
}

# --- 5. Pairwise comparison using consensus calls ---
# For each species pair, use only HOGs where BOTH species have a clean A/B
# consensus (excluding "Divergent" and HOGs missing from either species).
# Note: this means each pair may draw on a slightly different HOG subset —
# that's fine here since the fix targets copy-pairing noise, not HOG-set
# consistency; add a core-HOG filter beforehand if you want both.

pairs <- combn(species_list, 2, simplify = FALSE)
comparison_list <- vector("list", length(pairs))

for (i in seq_along(pairs)) {
  s1 <- pairs[[i]][1]
  s2 <- pairs[[i]][2]
  
  sub <- wide[
    !is.na(wide[[s1]]) & !is.na(wide[[s2]]) &
      wide[[s1]] != "Divergent" & wide[[s2]] != "Divergent",
    c("HOG", s1, s2)
  ]
  names(sub) <- c("HOG", "compartment1", "compartment2")
  sub$spc1 <- s1
  sub$spc2 <- s2
  sub$compartment_switch <- sub$compartment1 != sub$compartment2
  
  comparison_list[[i]] <- sub
}

pairwise_comparison <- do.call(rbind, comparison_list)

# --- 6. Summarise: % conserved vs. % switched per species pair ---

pairwise_summary <- aggregate(
  compartment_switch ~ spc1 + spc2,
  data = pairwise_comparison,
  FUN = function(x) c(n = length(x), n_switch = sum(x))
)

pairwise_summary <- do.call(data.frame, pairwise_summary)
names(pairwise_summary) <- c("spc1", "spc2", "n_comparisons", "n_switch")
pairwise_summary$n_conserved   <- pairwise_summary$n_comparisons - pairwise_summary$n_switch
pairwise_summary$pct_conserved <- pairwise_summary$n_conserved / pairwise_summary$n_comparisons
pairwise_summary$pct_switch    <- pairwise_summary$n_switch    / pairwise_summary$n_comparisons

# --- 7. Mirror to a full symmetric matrix ---

mirrored <- pairwise_summary
mirrored$spc1 <- pairwise_summary$spc2
mirrored$spc2 <- pairwise_summary$spc1

matrix_data <- rbind(pairwise_summary, mirrored)

# --- 8. Heatmap ---

library(ggplot2)

species_order <- c("Eluc", "Upyg", "Omyk", "Tthy", "Salp", "Ssal")
matrix_data$spc1 <- factor(matrix_data$spc1, levels = species_order)
matrix_data$spc2 <- factor(matrix_data$spc2, levels = rev(species_order))

conservation_heatmap <- ggplot(matrix_data, aes(x = spc1, y = spc2, fill = pct_conserved)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = scales::percent(pct_conserved, accuracy = 1)),
            color = "black", size = 3.5) +
  scale_fill_gradient(low = "#fee0d2", high = "#08519c",
                      limits = c(0, 1), name = "% compartment\nconserved") +
  labs(x = NULL, y = NULL,
       title = "Pairwise compartment conservation (consensus calls, no copy cross-join)") +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(face = "bold", size = 13)
  ) +
  coord_fixed()

print(pairwise_summary)
print(conservation_heatmap)

library(ape)

# --- 1. Build a square distance matrix from pairwise_summary ---

species_list <- c("Ssal", "Salp", "Omyk", "Tthy", "Eluc", "Upyg")
n <- length(species_list)

dist_matrix <- matrix(0, nrow = n, ncol = n,
                      dimnames = list(species_list, species_list))

for (i in seq_len(nrow(pairwise_summary))) {
  s1 <- pairwise_summary$spc1[i]
  s2 <- pairwise_summary$spc2[i]
  d  <- pairwise_summary$pct_switch[i]
  dist_matrix[s1, s2] <- d
  dist_matrix[s2, s1] <- d   # symmetric
}

# diagonal (species vs itself) stays 0, which is correct for a distance matrix

# --- 2. Convert to a formal "dist" object ---

compartment_dist <- as.dist(dist_matrix)

# --- 3. Build the tree ---
# Two common choices:
#   - UPGMA (hclust, method = "average"): assumes a constant rate of change,
#     produces an ultrametric (equal-depth) tree — simpler, good first look
#   - Neighbor-joining (ape::nj): does not assume constant rate, more
#     standard for "phylogeny-style" trees, allows unequal branch lengths

# UPGMA
upgma_tree <- hclust(compartment_dist, method = "average")

# Neighbor-joining
nj_tree <- nj(compartment_dist)

# --- 4. Plot ---

# UPGMA dendrogram (base R)
plot(upgma_tree, main = "Compartment divergence, consensus calls (UPGMA)",
     xlab = "", sub = "", ylab = "Switch rate distance")

# Neighbor-joining tree (ape)
plot(nj_tree, main = "Compartment divergence, consensus calls (neighbor-joining)", type = "unrooted")
add.scale.bar()

# bootstraping
species_list <- c("Ssal", "Salp", "Omyk", "Tthy", "Eluc", "Upyg")

# --- 1. Restrict to "complete case" HOGs: every species has a clean A/B call ---
# (no NA, no "Divergent") — this is the character matrix the bootstrap resamples

is_clean <- apply(wide[, species_list], 1, function(row) all(row %in% c("A", "B")))
complete_wide <- wide[is_clean, ]

cat(sprintf("Complete-case HOGs used for bootstrap: %d (out of %d total HOGs)\n",
            nrow(complete_wide), nrow(wide)))

# character matrix: rows = species (taxa), columns = HOGs (characters) —
# same orientation ape::boot.phylo expects for sequence alignments
char_matrix <- t(as.matrix(complete_wide[, species_list]))
rownames(char_matrix) <- species_list

# --- 2. Function: build an NJ tree from a (possibly resampled) character matrix ---

build_switch_tree <- function(x) {
  sp <- rownames(x)
  n  <- length(sp)
  d  <- matrix(0, n, n, dimnames = list(sp, sp))
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      switch_rate <- mean(x[i, ] != x[j, ])
      d[i, j] <- switch_rate
      d[j, i] <- switch_rate
    }
  }
  nj(as.dist(d))
}

# --- 3. Build the tree on the real (non-resampled) data ---

original_tree <- build_switch_tree(char_matrix)

# --- 4. Bootstrap: resample HOGs with replacement, rebuild tree, repeat ---

set.seed(1)
bootstrap_support <- boot.phylo(
  original_tree, char_matrix, build_switch_tree,
  B = 1000, rooted = FALSE
)

# --- 5. Plot with bootstrap support values at each node ---

plot(original_tree, type = "unrooted",
     main = "Compartment divergence, NJ tree with bootstrap support")
nodelabels(bootstrap_support, frame = "none", adj = c(1.2, -0.3), cex = 0.9)
add.scale.bar()
