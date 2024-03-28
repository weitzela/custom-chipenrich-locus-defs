##### Prepare Environment #####
pacman::p_load(GenomicRanges, tidyverse)

if (basename(getwd()) == "R") {
  setwd("..")
}

# define genome specific objects
genome_id = "rn7"
ensembl = biomaRt::useEnsembl(biomart = "genes", dataset = "rnorvegicus_gene_ensembl", version = 109)

# Obtain all gene IDs that exist in the pathway databases probed, to get as many 
# ensembl ID matches as possible. These entrez IDs are used to match back to ensembl 
# IDs through biomaRt because it does not return consistent matches from ensembl 
# IDs to entrez, and vice versa.
all_entrez_ids = c(chipenrich.data::geneset.GOBP.rno@all.genes,
                   chipenrich.data::geneset.GOCC.rno@all.genes,
                   chipenrich.data::geneset.GOMF.rno@all.genes,
                   chipenrich.data::geneset.kegg_pathway.rno@all.genes) %>%
  unique()

# load functions
source("R/fxns.R")

##### Create Custom Locus Definitions #####
anno_gr = prepareTSSandGeneBodyGR()
anno_gr$promoter = GenomicRanges::promoters(anno_gr$tss, upstream = 2000, downstream = 200)

locus_defs = keep_at(anno_gr, c("promoter", "gene_body"))
locus_defs$promoter = makeLocusDefsAroundTSS(anno_gr$tss, anno_gr$promoter)
locus_defs$intergenic = makeBpOutsideLocusDef(anno_gr$tss, sort(unlist(locus_defs, use.names = FALSE)))
locus_defs[["2kb_outside"]] = makeBpOutsideLocusDef(anno_gr$tss, 2000)
locus_defs[["2kb_around_TSS"]] = makeLocusDefsAroundTSS(anno_gr$tss, 2000)

# get ensembl:gene ID matches to convert ensembl IDs to entrez IDs
ensembl_entrez_pairs = list(
  "ensembl_gene_id" = anno_gr$tss$ensembl_gene_id, 
  "entrezgene_id" = all_entrez_ids
) %>%
  imap(~ biomaRt::getBM(attributes = c("entrezgene_id", "ensembl_gene_id"), filters = .y, values = .x, mart = ensembl)) %>%
  bind_rows() %>%
  distinct() %>%
  drop_na() %>%
  dplyr::rename("gene_id" = "entrezgene_id")

locus_defs = imap(locus_defs, function(.x, .y) {
  .x %>%
    `seqlevelsStyle<-`("UCSC") %>%
    as.data.frame() %>%
    select(1:3, ensembl_gene_id) %>% 
    left_join(ensembl_entrez_pairs, ., by = "ensembl_gene_id", relationship = "many-to-many") %>%
    arrange(seqnames, start) %>%
    relocate(gene_id, .after = last_col()) %>%
    drop_na(seqnames) %>%
    select(-ensembl_gene_id) %>%
    `attr<-`("filename", paste0(genome_id, "-locus_", str_replace_all(.y, "_", "-"), ".txt"))
})

walk(locus_defs, ~ write_tsv(.x, file.path("locus-defs", attr(.x, "filename"))))