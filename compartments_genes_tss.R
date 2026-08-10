library(tidyverse)
library(rtracklayer)
library(GenomicRanges)

genes_gff <- import(
  "genomes/Salmo_salar.Ssal_v3.1.115.chr.gff3",
  feature.type = "gene"
)
tss_ssal <- resize(genes_gff, width = 1, fix = "start")

tss_df <- tibble(
  chrom = as.character(seqnames(tss_ssal)),
  tss_start = start(tss_ssal),
  tss_end = end(tss_ssal),
  strand = as.character(strand(tss_ssal)),
  gene_id = mcols(tss_ssal)$gene_id
)

comp_liver_ssal <- read_tsv("compartment_liver_50kb/ab_compartment_liver_gene_density_Ss.50kb.tsv", show_col_types = FALSE)

comp_gr <- GRanges(
  seqnames = comp_liver_ssal$chrom,
  ranges = IRanges(start = comp_liver_ssal$start + 1, end = comp_liver_ssal$end),
  E1 = comp_liver_ssal$E1,
  E2 = comp_liver_ssal$E2,
  E3 = comp_liver_ssal$E3,
  compartment = comp_liver_ssal$compartment
)

hits <- findOverlaps(tss_ssal, comp_gr)

tss_with_comp <- tss_df[queryHits(hits), ] %>%
  mutate(
    E1 = mcols(comp_gr)$E1[subjectHits(hits)],
    E2 = mcols(comp_gr)$E2[subjectHits(hits)],
    E3 = mcols(comp_gr)$E3[subjectHits(hits)],
    compartment = mcols(comp_gr)$compartment[subjectHits(hits)]
  )

# tss is only 1-bp genomic interval, single genomic coordinate, the first bp of the gene

ssal_genes <- read_tsv("allDup_withinSpecies_busco.tsv", show_col_types = FALSE) |> 
  filter(spc == "Ssal")


gene_assigned_comp <- ssal_genes |> 
  left_join(
    tss_with_comp |> 
      select(gene_id, chrom, tss_start, tss_end, strand, E1, E2, E3, compartment), 
    by = "gene_id"
  )

write_tsv(gene_assigned_comp, file = "genes_assigned_compartment/Ssal_genes_with_TSS_compartments.tsv")

