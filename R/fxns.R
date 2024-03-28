assignMidpointRegionsToGenes = function(tss_gr, gene_id_col = "ensembl_gene_id") {
  # make sure tss strand is *
  strand(tss_gr) = "*"
  # get gaps between input TSSs
  gaps_btwn_tss = GenomicRanges::gaps(tss_gr) %>% .[(strand(.) == "*")]
  # get midpoints of those gaps, which will serve as boundaries between genes
  btwn_gene_barriers = gUtils::gr.mid(gaps_btwn_tss)
  # fill in the space between those barriers to begin definition of locus defs
  locus_def = GenomicRanges::gaps(btwn_gene_barriers) %>% .[(strand(.) == "*")]
  # match the tss that overlaps the locus definition
  tss_locus_overlap = GenomicRanges::findOverlaps(locus_def, tss_gr)
  locus_def = locus_def[tss_locus_overlap@from]
  mcols(locus_def)[[gene_id_col]] = mcols(tss_gr[tss_locus_overlap@to])[[gene_id_col]]
  return(locus_def)
}

makeBpOutsideLocusDef = function(tss_gr, region_to_rm = NULL, gene_id_col = "ensembl_gene_id") {
  if (is.numeric(region_to_rm)) {
    region_to_rm = flank(tss_gr, region_to_rm, both = TRUE)
  } else if ((!inherits(region_to_rm, "GRanges")) | is.null(region_to_rm)) {
    stop("Must define variable 'region_to_rm' as either a numberic value indicating distance around TSS or GRanges object.")
  }
  locus_def = assignMidpointRegionsToGenes(tss_gr, gene_id_col)
  # remove X around TSS for each gene
  locus_def = GenomicRanges::subtract(locus_def, region_to_rm, ignore.strand = TRUE) %>% 
    `names<-`(mcols(locus_def)[[gene_id_col]]) %>% 
    unlist()
  mcols(locus_def)[[gene_id_col]] = names(locus_def)
  locus_def = sort(locus_def) %>% 
    `names<-`(NULL)
  return(locus_def)
}

makeLocusDefsAroundTSS = function(tss_gr, region_to_keep = NULL, gene_id_col = "ensembl_gene_id") {
  if (is.numeric(region_to_keep)) {
    region_to_keep = flank(tss_gr, region_to_keep, both = TRUE)
  } else if ((!inherits(region_to_keep, "GRanges")) | is.null(region_to_keep)) {
    stop("Must define variable 'region_to_keep' as either a numberic value indicating distance around TSS or GRanges object.")
  }
  join_cols = map_chr(c(".x", ".y"), ~ paste0(gene_id_col, .x))
  gr = assignMidpointRegionsToGenes(tss_gr, gene_id_col) %>% 
    plyranges::join_overlap_intersect(., region_to_keep) %>% 
    plyranges::filter(!!rlang::sym(join_cols[1]) == !!rlang::sym(join_cols[2])) %>% 
    plyranges::mutate("{gene_id_col}" := !!rlang::sym(join_cols[1])) %>% 
    plyranges::select(-any_of(join_cols))
  return(gr)
}

prepareTSSandGeneBodyGR = function() {
  genome_info = readRDS(file.path("genome-info", paste0(genome_id, "_seq-info.RData")))
  all_ensembls = biomaRt::getBM(attributes = "ensembl_gene_id", mart = ensembl) %>%
    pull(ensembl_gene_id)
  gr_ls = list("tss" = c("transcription_start_site", "transcript_is_canonical"),
               "gene_body" = c("start_position", "end_position")) %>%
    imap(function(.biomart_info, .anno_label) {
      biomart_out = biomaRt::getBM(attributes = c(.biomart_info, "ensembl_gene_id", "chromosome_name", "strand"), 
                                   values = all_ensembls,
                                   mart = ensembl) %>%
        mutate(strand = ifelse(strand == 1, "+", "-")) %>%
        na.omit()
      if (.anno_label == "tss") {
        biomart_out = biomart_out %>%
          dplyr::rename("tss" = "transcription_start_site") %>%
          mutate(start = tss, end = tss) %>%
          select(-transcript_is_canonical)
      } else if (.anno_label == "gene_body") {
        biomart_out = biomart_out %>%
          dplyr::rename_with(function(.c) str_remove_all(.c, "_position"))
      }
      gr = makeGRangesFromDataFrame(biomart_out, keep.extra.columns = TRUE) %>% 
        subset(seqnames(.) %in% c(1:20))
      seqlevels(gr) = seqlevelsInUse(gr)
      genome(gr) = genome(genome_info)
      GenomicRanges::seqinfo(gr) = genome_info
      gr = sort(gr)
      names(gr) = NULL
      return(gr)
    }) %>%
    GRangesList()
  return(gr_ls)
}

