library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
# looped - genes assigned liver compartments

# Read once, then filter inside the loop.
all_dup <- read_tsv("allDup_withinSpecies_busco.tsv", show_col_types = FALSE)

samples <- tribble(
  ~spc,    ~gff_file,                                                   ~comp_file,
  "Eluc",  "genomes/Esox_lucius.fEsoLuc1.pri.115.chr.gff3",             "compartment_liver_50kb/ab_compartment_L11_gene_density_El.50kb.tsv",
  "Salp",  "genomes/Salvelinus_alpinus_Arthur.gff3",                    "compartment_liver_50kb/ab_compartment_liver_gene_density_Sa.50kb.tsv",
  "Ssal",  "genomes/Salmo_salar.Ssal_v3.1.115.chr.gff3",                "compartment_liver_50kb/ab_compartment_liver_gene_density_Ss.50kb.tsv",
  "Tthy", "genomes/Thymallus_thymallus.Garry.chr.gff3",                 "compartment_liver_50kb/ab_compartment_liver_gene_density_Tt.50kb.tsv",
  "Upyg", "genomes/GCA_040894045.1_UPYG_1.0_genomic.chr.gff",           "compartment_liver_50kb/ab_compartment_liver_gene_density_Up.50kb.tsv",
  "Omyk",  "genomes/Oncorhynchus_mykiss.USDA_OmykA_1.1.115.chr.gff3",   "compartment_liver_50kb/ab_compartment_liver_gene_density_Om.50kb.tsv"
)

assign_compartments <- function(spc, gff_file, comp_file, all_dup) {
  message("Processing:", spc)
  
  genes_gff <- import(gff_file, feature.type = "gene")
  tss_gr <- resize(genes_gff, width =1, fix = "start")
  
  gene_metadata <- mcols(tss_gr)
  
  if ("gene_id" %in% colnames(gene_metadata)) {
    gff_gene_id <- as.character(gene_metadata$gene_id)
    
  } else if ("ID" %in% colnames(gene_metadata)) {
    gff_gene_id <- as.character(gene_metadata$ID)
    
  } else {
    stop("Neither gene_id nor ID was found in: ", gff_file)
  }
  
  tss_df <- tibble(
    chrom = as.character(seqnames(tss_gr)),
    tss_start = start(tss_gr),
    tss_end = end(tss_gr),
    strand = as.character(strand(tss_gr)),
    gene_id = gff_gene_id
  )
  
  comp <- read_tsv(comp_file, show_col_types = FALSE)
  
  comp_gr <- GRanges(
    seqnames = comp$chrom,
    ranges = IRanges(start = comp$start + 1, end = comp$end),
    E1 = comp$E1,
    E2 = comp$E2,
    E3 = comp$E3,
    compartment = comp$compartment
  )
  
  hits <- findOverlaps(tss_gr, comp_gr)
  
  tss_with_comp <- tss_df[queryHits(hits), ] %>%
    mutate(
      E1 = mcols(comp_gr)$E1[subjectHits(hits)],
      E2 = mcols(comp_gr)$E2[subjectHits(hits)],
      E3 = mcols(comp_gr)$E3[subjectHits(hits)],
      compartment = mcols(comp_gr)$compartment[subjectHits(hits)]
    )
  
  gene_assigned_comp <- all_dup |> 
    filter(spc == !!spc) |> 
    left_join(
      tss_with_comp |> 
        select(gene_id, chrom, tss_start, tss_end, strand, E1, E2, E3, compartment), 
      by = "gene_id"
    )
  
  output_file <- file.path("genes_assigned_compartment", paste0(spc, "_genes_wTSS_compartment_liver.tsv"))
  
  write_tsv(gene_assigned_comp, output_file)
  
  message(
    spc, ":", sum(!is.na(gene_assigned_comp$compartment)),
    "/", nrow(gene_assigned_comp),
    "allDup genes received a compartment call."
  )
  
  gene_assigned_comp
  
}

results <- purrr::pmap(
  samples, 
  ~ assign_compartments(
    spc = ..1,
    gff_file = ..2,
    comp_file = ..3,
    all_dup = all_dup
  )
)
