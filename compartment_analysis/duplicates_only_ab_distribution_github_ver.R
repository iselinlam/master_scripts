library(tidyverse)
library(dplyr)

spc_levels    <- c("Tthy", "Salp", "Ssal", "Omyk")
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

combos <- expand_grid(spc = spc_levels, tissue = tissue_levels)

all_genes <- read_tsv("all_genes_SS4R.tsv")

# -- all_genes_SS4R.tsv filtered with this code from a total list of duplicate genes with compartments assigned by the TSS.
# ss4ronly <- read_tsv(file[1], show_col_types = F) |> 
# mutate(spc = spc_i, tissue = tissue_i) |> 
#   filter(Type == "SS4Ronly", compartment %in% c("A", "B")) |> 
#   add_count(Orthogroup, name = "n_copies") |>
#   filter(n_copies == 2) |>
#   select(-n_copies, -seqname, -BeforeSS4R, -AfterSS4R, -Outgroup, -Tip, -tss_start, -tss_end)


species_per_tissue <- all_genes |> 
  distinct(tissue, spc) |> 
  count(tissue, name = "n_expected_species")

ss4ronly_genes <- all_genes |> 
  inner_join(species_per_tissue, by = "tissue") |> 
  group_by(tissue, Orthogroup) |>
  filter(
    n_distinct(spc) == dplyr::first(n_expected_species),
    n() == 2 * dplyr::first(n_expected_species)
  ) |>
  ungroup() |> 
  select(-n_expected_species)

ss4r_count <- ss4ronly_genes |> 
  group_by(spc, tissue, Orthogroup) |> 
  summarise(
    copy1 = compartment[1],
    copy2 = compartment[2],
    .groups = "drop"
  ) %>%
  mutate(
    pair_state = case_when(
      copy1 == "A" & copy2 == "A" ~ "AA",
      copy1 == "B" & copy2 == "B" ~ "BB",
      TRUE ~ "AB"
    )
  )

ss4r_ab_distribution <- ss4r_count |> 
  group_by(spc, tissue) |> 
  count(pair_state) |> 
  mutate(pct = round( n /sum(n), 3), spc = spc, tissue = tissue) |>
  select(-n) |> 
  pivot_wider(
    names_from = pair_state, 
    values_from = pct, 
    names_prefix = "obs_"
  )

ab_genomewide <- read_tsv("ab_fractions_genes_genomewide.tsv") |> 
  filter(spc %in% spc_levels) 

ab_duplicates_distribution <- ab_genomewide |> 
  left_join(
    ss4r_ab_distribution, 
    by = c("spc", "tissue")
  )

plot_data <- ab_duplicates_distribution |>
  select(spc, tissue, matches("^(exp|obs)_(AA|AB|BB)$")) |>
  pivot_longer(
    cols = matches("^(exp|obs)_(AA|AB|BB)$"),
    names_to = c("source", "class"),
    names_pattern = "^(exp|obs)_(AA|AB|BB)$",
    values_to = "proportion"
  ) |>
  mutate(
    source = recode(source, exp = "Expected", obs = "Observed"),
    class = factor(class, levels = c("AA", "AB", "BB")),
  ) |>
  group_by(spc, tissue, source) |>
  mutate(proportion = proportion / sum(proportion, na.rm = TRUE)) |>
  ungroup()

plot_data_for_figures <- plot_data |>
  filter(!(spc == "Tthy" & tissue == "gill")) |> 
  left_join(
    species_names |> select(spc, common_name), 
    by = "spc"
  )

plots_by_spc <- split(
  plot_data_for_figures,
  plot_data_for_figures$spc
) |>
  lapply(function(dat) {
    
    dat <- dat |>
      mutate(
        source = factor(source,
                        levels = c("Observed", "Expected"),
                        labels = c("Obs", "Exp")),
        tissue = factor(tissue, levels = c("brain", "gill", "liver"))
      )
    
    ggplot(dat, aes(x = source, y = proportion, fill = class)) +
      geom_col(
        width = 0.78,
        colour = "black",
        linewidth = 0.3,
        position = position_stack(reverse = TRUE)
      ) +
      geom_text(
        aes(label = scales::percent(proportion, accuracy = 0.1)),
        position = position_stack(vjust = 0.5, reverse = TRUE),
        colour = "white",
        fontface = "bold",
        size = 4
      ) +
      facet_grid(. ~ tissue, scales = "free_x", space = "free_x") +
      scale_fill_manual(
        values = c(AA = "#B52A00", AB = "#7656A6", BB = "#2F6F9F"),
        name = "Pair type"
      ) +
      scale_y_continuous(
        labels = scales::percent,
        limits = c(0, 1),
        expand = expansion(mult = c(0, 0.03))
      ) +
      labs(
        title = paste0(unique(dat$spc), " - ", unique(dat$common_name)),
        x = NULL,
        y = "Fraction of gene pairs (%)"
      ) +
      theme_classic(base_size = 14) +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 14),
        legend.position = "bottom"
      )
  })

purrr::iwalk(
  plots_by_spc,
  ~ ggsave(
    path = "distribution",
    filename = paste0(.y, "_duplicates_distribution.png"),
    plot = .x,
    width = 9,
    height = 6,
    dpi = 300
  )
)


