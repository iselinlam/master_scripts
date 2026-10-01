library(tidyverse)

spc_levels    <- c("Tthy", "Salp", "Ssal", "Omyk", "Eluc", "Upyg")
tissue_levels <- c("liver", "brain", "gill")

species_names <- tribble(
  ~spc, ~common_name, ~scientific, 
  "Ssal", "Atlantic Salmon", "Salmo salar",
  "Salp", "Artic Char", "Salvelinus_alpinus",
  "Omyk", "Rainbow Trout", "Oncorhynchus mykiss",
  "Tthy", "European Grayling", "Thymallus thymallus", 
  "Eluc", "Northern Pike", "Esox lucius",
  "Upyg", "Eastern Mudminnow", "Umbra pygmaea"
)

# ---- input file made with this code: 

# combos <- expand_grid(spc = spc_levels, tissue = tissue_levels)
# 
# files_list <- list()
# 
# for (i in seq_len(nrow(combos))) {
#   spc_i    <- combos$spc[i]
#   tissue_i <- combos$tissue[i]
#   file <- file.path("all_duplicates_with_compartments",
#                     paste0(spc_i, "_genes_wTSS_compartment_", tissue_i, ".tsv"))
#   
#   if (!file.exists(file)) {
#     message(spc_i, " - ", tissue_i, ": file not found, skipping")
#     next
#   }
#   
#   singleton_genes <- read_tsv(file, show_col_types = FALSE) |> 
#     mutate(spc = spc_i, tissue = tissue_i) |> 
#     filter(Type == "Singleton") |> 
#     add_count(Orthogroup, name = "n_copies") |>   # copies within this species
#     filter(n_copies == 1) |>
#     filter(compartment %in% c("A", "B")) |>       # compartment filter after copy count
#     select(-n_copies, -seqname, -BeforeSS4R, -AfterSS4R, -Outgroup, -Tip,
#            -tss_start, -tss_end)
#   
#   files_list[[paste(spc_i, tissue_i, sep = "_")]] <- singleton_genes
# }
# 
# singleton_combined <- bind_rows(
#   map(files_list, ~mutate(.x, chrom = as.character(chrom))))
# 
# write_tsv(singleton_combined, "singleton_analysis/singleton_combined.tsv", show_col_types = FALSE)


# ---- input file:
singleton_combines <- read_tsv("singleton_combined.tsv)

species_per_tissue <- singleton_combined |>
  distinct(tissue, spc) |>
  group_by(tissue) |>
  summarise(n_spc = n(), spcs = paste(sort(spc), collapse = ", "))
print(species_per_tissue)

singleton_1to1 <- singleton_combined |>
  left_join(species_per_tissue |> select(tissue, n_spc), by = "tissue") |>
  group_by(tissue, Orthogroup) |>
  filter(n_distinct(spc) == dplyr::first(n_spc)) |>
  ungroup() |>
  select(-n_spc)

singleton_1to1 |>
  distinct(tissue, Orthogroup) |>
  count(tissue, name = "n_orthogroups") |>
  print()

walk(tissue_levels, function(t) {
  singleton_1to1 |>
    filter(tissue == t) |>
    write_tsv(file.path("singleton_analysis",
                        paste0("singleton_1to1_", t, ".tsv")))
})

## pairwise and PLOT

species_order_all <- c("Upyg", "Eluc", "Tthy", "Omyk", "Salp", "Ssal")
run_tissue_heatmap <- function(tissue) {
  
  # --- Read 1:1 singleton file and pivot to wide A/B calls ---
  dat <- read_tsv(file.path("singleton_analysis",
                            paste0("singleton_1to1_", tissue, ".tsv")),
                  show_col_types = FALSE)
  
  compartment_wide <- dat |>
    select(Orthogroup, spc, compartment) |>
    distinct() |>
    pivot_wider(names_from = spc, values_from = compartment)
  
  species_list  <- setdiff(colnames(compartment_wide), "Orthogroup")
  species_order <- intersect(species_order_all, species_list)
  species_pairs <- combn(species_list, 2, simplify = FALSE)
  
  # --- Pairwise % compartment identity ---
  pairwise_conservation <- map_dfr(species_pairs, function(pair) {
    sp1 <- pair[1]
    sp2 <- pair[2]
    
    pair_data <- compartment_wide |>
      select(Orthogroup, all_of(sp1), all_of(sp2)) |>
      drop_na()
    
    n_pairs <- nrow(pair_data)
    
    if (n_pairs == 0) {
      return(tibble(species1 = sp1, species2 = sp2, n_orthologs = 0L,
                    n_conserved = NA_integer_, pct_conserved = NA_real_,
                    p_value = NA_real_))
    }
    
    n_conserved <- sum(pair_data[[sp1]] == pair_data[[sp2]])
    binom_test  <- binom.test(n_conserved, n_pairs, p = 0.5, alternative = "greater")
    
    tibble(species1 = sp1, species2 = sp2,
           n_orthologs = n_pairs, n_conserved = n_conserved,
           pct_conserved = 100 * n_conserved / n_pairs,
           p_value = binom_test$p.value)
  }) |>
    mutate(p_adj = p.adjust(p_value, method = "BH"),
           tissue = tissue)
  
  write_tsv(pairwise_conservation,
            file.path("singleton_analysis",
                      paste0("singleton_compartment_conservation_", tissue, ".tsv")))
  
  # --- Heatmap ---
  heatmap_data <- pairwise_conservation |>
    select(species1, species2, pct_conserved, n_orthologs) |>
    bind_rows(
      pairwise_conservation |>
        select(species1 = species2, species2 = species1, pct_conserved, n_orthologs)
    ) |>
    mutate(species1 = factor(species1, levels = species_order),
           species2 = factor(species2, levels = species_order))
  
  p <- ggplot(heatmap_data, aes(x = species2, y = species1, fill = pct_conserved)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.1f%%", pct_conserved)),
              color = "black", size = 4.2, fontface = "bold", vjust = -0.3) +
    geom_text(aes(label = paste0("n=", n_orthologs)),
              color = "black", size = 2.8, alpha = 0.85, vjust = 1.6) +
    scale_fill_gradient2(low = "#3B4CC0", mid = "#DDDDDD", high = "#B40426",
                         midpoint = 50, limits = c(0, 100), name = "% conserved") +
    labs(title = paste("Singleton ortholog A/B compartment conservation -", tissue),
         subtitle = "Pairwise % identical compartment call") +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(), axis.title = element_blank())
  
  ggsave(file.path("singleton_analysis",
                   paste0("singleton_compartment_conservation_", tissue, ".png")),
         p, width = 8, height = 6, dpi = 200)
  
  list(data = pairwise_conservation, plot = p)
}

tissue_levels <- c("liver", "brain", "gill")
results <- set_names(tissue_levels) |> map(run_tissue_heatmap)

results$brain$plot   # view any tissue's plot


library(ape)

outgroup_all <- c("Eluc", "Upyg")
name_map <- setNames(species_names$common_name, species_names$spc)

run_tissue_trees <- function(tissue, B = 1000) {
  
  # --- 1. Distance matrix from pairwise % conservation ---
  pw <- read_tsv(file.path("singleton_analysis",
                           paste0("singleton_compartment_conservation_", tissue, ".tsv")),
                 show_col_types = FALSE)
  
  species_list <- sort(unique(c(pw$species1, pw$species2)))
  n <- length(species_list)
  dist_matrix <- matrix(0, n, n, dimnames = list(species_list, species_list))
  
  for (i in seq_len(nrow(pw))) {
    d <- 1 - pw$pct_conserved[i] / 100
    dist_matrix[pw$species1[i], pw$species2[i]] <- d
    dist_matrix[pw$species2[i], pw$species1[i]] <- d
  }
  compartment_dist <- as.dist(dist_matrix)
  
  # --- 2. UPGMA and NJ trees ---
  upgma_tree <- hclust(compartment_dist, method = "average")
  nj_tree    <- nj(compartment_dist)
  
  # Root on outgroups only if they are monophyletic in the NJ tree
  outgroup    <- intersect(outgroup_all, species_list)
  use_rooting <- length(outgroup) == 2 && is.monophyletic(nj_tree, outgroup)
  message(tissue, ": outgroups monophyletic in NJ tree = ", use_rooting,
          " -> ", if (use_rooting) "rooting on outgroups" else "unrooted")
  
  nj_final <- if (use_rooting) root(nj_tree, outgroup = outgroup, resolve.root = TRUE) else nj_tree
  
  # --- 3. Per-ortholog compartment matrix for bootstrapping ---
  dat <- read_tsv(file.path("singleton_analysis",
                            paste0("singleton_1to1_", tissue, ".tsv")),
                  show_col_types = FALSE)
  compartment_wide <- dat |>
    select(Orthogroup, spc, compartment) |>
    distinct() |>
    pivot_wider(names_from = spc, values_from = compartment)
  
  compartment_matrix <- t(as.matrix(compartment_wide[, species_list]))
  compartment_matrix <- compartment_matrix[, colSums(is.na(compartment_matrix)) == 0]
  cat(sprintf("%s: %d orthologs used for bootstrap\n", tissue, ncol(compartment_matrix)))
  
  # --- 4. Bootstrap ---
  build_match_tree <- function(x) {
    sp <- rownames(x); ns <- length(sp)
    d <- matrix(0, ns, ns, dimnames = list(sp, sp))
    for (i in seq_len(ns - 1)) for (j in (i + 1):ns) {
      v <- 1 - mean(x[i, ] == x[j, ])
      d[i, j] <- v; d[j, i] <- v
    }
    nj(as.dist(d))
  }
  
  set.seed(1)
  boot_result <- boot.phylo(nj_tree, compartment_matrix, build_match_tree,
                            B = B, rooted = FALSE, trees = TRUE)
  
  if (use_rooting) {
    boot_rooted <- lapply(boot_result$trees, function(t) {
      tryCatch(root(t, outgroup = outgroup, resolve.root = TRUE),
               error = function(e) NULL)
    })
    boot_rooted <- Filter(Negate(is.null), boot_rooted)
    class(boot_rooted) <- "multiPhylo"
    support <- prop.clades(nj_final, boot_rooted, rooted = TRUE)  # counts out of B
    cat(sprintf("%s: %d/%d replicates had monophyletic outgroups\n",
                tissue, length(boot_rooted), B))
  } else {
    support <- prop.clades(nj_tree, boot_result$trees, rooted = FALSE)  # counts out of B
  }
  support[is.na(support)] <- 0
  
  # --- 5. Relabel tips with common names and plot ---
  nj_plot <- nj_final
  nj_plot$tip.label <- unname(name_map[nj_plot$tip.label])
  
  upgma_plot <- upgma_tree
  upgma_plot$labels <- unname(name_map[upgma_plot$labels])
  
  graphics.off()
  pdf(file.path("singleton_analysis",
                paste0("compartment_nj_upgma_trees_", tissue, ".pdf")),
      width = 12, height = 5)
  par(mfrow = c(1, 2), mar = c(4, 2, 3, 6))
  
  if (use_rooting) {
    plot(nj_plot, main = paste0("Rooted NJ (outgroup: Pike + Mudminnow) - ", tissue),
         cex = 0.9, label.offset = 0.01, y.lim = c(-1, Ntip(nj_plot)), font = 2)
    nodelabels(support, frame = "none", adj = c(0.5, -0.8), cex = 0.8, font = 2)
    add.scale.bar(x = 0, y = -0.5, cex = 0.8)
  } else {
    plot(nj_plot, type = "unrooted", main = paste0("Unrooted NJ - ", tissue),
         cex = 0.9, font = 2, no.margin = FALSE)
    nodelabels(support, frame = "none", cex = 0.8, font = 2)
    add.scale.bar(cex = 0.8)
  }
  mtext(sprintf("Node labels: bootstrap support (out of %d replicates)", B),
        side = 1, line = 2.5, cex = 0.8)
  
  plot(upgma_plot, main = paste0("UPGMA (A/B compartment) - ", tissue),
       xlab = "", sub = "", ylab = "1 - % conserved", cex = 0.9)
  
  par(mfrow = c(1, 1))
  dev.off()
  
  invisible(list(dist = compartment_dist, nj = nj_final, upgma = upgma_tree,
                 support = support, rooted = use_rooting))
}

tissue_levels <- c("liver", "brain", "gill")
tree_results <- set_names(tissue_levels) |> map(run_tissue_trees)




