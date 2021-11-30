#' Check gene list
#'
#' @param gff gff file
#' @param quename name of focal species
#' @param gene_tag a list of genes as character vector
#' @export
#' @return warning message
gene_check <- function(gff, quename, gene_tag) {

  "Assert query and control ids match the names in the GFF file. Also confirm
  that query and control ids do not overlap"

  check_ids <- function(ids){
    txnames <- GenomicRanges::mcols(gff)$attr
    matches <- ids %in% txnames
    if(!all(matches)){
      msg <- "%s of %s ids are missing in the focal GFF (%s). Here is the head: [%s]"
      stop(sprintf(
        msg,
        sum(!matches),
        length(ids),
        quename,
        paste(head(ids[!matches]), collapse=", ")
      ))
    }
  }

  check_ids(gene_tag)
}


#' Do alignment for the genes from focal species agaginst to target genome
#'
#' @param con configuration
#' @param quesp a list of species metadata for focal species from load_species
#' @param tarsp a list of species metadata for target species from load_species
#' @param quename a string of focal species name
#' @param tarname a string of target species name
#' @param gene_tag a list of genes as character vector
#' @export
#' @return a list of alignment result
compare_target_to_focal <- function(con, quesp, tarsp, quename, tarname, gene_tag){

  message("Processing ", quename, " vs ", tarname)

  qseqinfo <- quesp$seqinfo
  tseqinfo <- tarsp$seqinfo

  synmap <- synder::read_synmap(con@input@syn[[quename]][[tarname]],
                                seqinfo_a=qseqinfo,
                                seqinfo_b=tseqinfo)

  synmap_summary <- summarize_syn(synmap)

  nsims_prot <- con@alignment@simulation@prot2prot
  nsims_allorf <- con@alignment@simulation@prot2allorf
  nsims_transorf <- con@alignment@simulation@prot2transorf
  subMat <- con@alignment@substitution_matrix

  # SyntenyMap -> mRNA -> SearchIntervals
  synder_out <- synder::search(
    syn     = synmap,
    gff     = quesp$mRNA,
    swap    = FALSE,
    trans   = con@synder@trans,
    k       = con@synder@k,
    r       = con@synder@r,
    tcl     = tseqinfo,
    qcl     = qseqinfo,
    offsets = con@synder@offsets
  )

  summary_synder_flags <- synder::flag_summary(synder_out)

  indels <- find_indels(synder_out, indel.threshold=con@alignment@indel_threshold)


  gaps <- overlapMap(gff=tarsp$nstring, si=synder_out) %>%
    {
      if(is.null(tarsp$nstring)){
        data.frame(
          query  = character(0),
          siid   = integer(0),
          length = integer(0)
        )
      } else {
        data.frame(
          query  = .$query,
          siid   = .$qid,
          length = GenomicRanges::width(tarsp$nstring)[.$qid]
        ) %>% { .[!is.na(.$length), ] }
      }
    }

  f_si_map <- overlapMap(gff=tarsp$mRNA, si=synder_out)

  f_si_map_orf <- overlapMap(gff=tarsp$orfgff, si=synder_out)

  aa2aa <- align_by_map(
    queseq = quesp$faa,
    tarseq = tarsp$faa,
    map = f_si_map,
    queries = gene_tag,
    nsims=nsims_prot,
    substitution_matrix=subMat
  )

  aa2orf <- align_by_map(
    queseq = quesp$faa,
    tarseq = tarsp$orffaa,
    map = f_si_map_orf,
    queries = gene_tag,
    nsims=nsims_allorf,
    substitution_matrix=subMat
  )

  f_si_map_transorf <- tarsp$transorfgff %>%
    {
      data.frame(
        seqid = GenomicRanges::seqnames(.),
        orfid = GenomicRanges::mcols(.)$attr,
        stringsAsFactors=FALSE
      )
    } %>%
    merge(f_si_map, by.x="seqid", by.y="target") %>%
    {
      data.frame(
        query = .$query,
        target = .$orfid,
        type = "transorf",
        stringsAsFactors=FALSE
      )
    }

  aa2transorf <- align_by_map(
    queseq = quesp$faa,
    tarseq = tarsp$transorfaa,
    map = f_si_map_transorf,
    queries = gene_tag,
    nsims=nsims_transorf,
    substitution_matrix=subMat
  )


  # gene2genome - align query DNA sequence against the SI

  map     = synder_out
  quedna  = quesp$transcriptomeSeq
  tardna  = tarsp$genomeDB
  queries = gene_tag

  qids <- GenomicRanges::mcols(map)$attr
  if(! all(qids %in% names(quedna)) ){

    ngood <- sum(unique(qids) %in% names(quedna))
    ntotal <- length(unique(qids))

    msg <- "There is a mismatch between the names in the synteny map and
      those derived from the GFF file. %s of %s names from the synteny map are
      missing in the DNA file. Here are the first 5 names from the GFF-derived
      DNA file: [%s]. Here are the first 5 from the synteny map: [%s]."

    stop(sprintf(
      msg,
      ngood,
      ntotal,
      paste0(head(names(quedna), 5), collapse=", "),
      paste0(head(qids, 5), collapse=", ")
    ))

  }
  
  queries <- queries[queries %in% qids]
  tarseq <- Rsamtools::getSeq(x=tardna, CNEr::second(map))
  queseq <- quedna[qids]
  offset <- GenomicRanges::start(CNEr::second(map)) - 1

  gene2genome_aln <- get_dna2dna(
    queseq   = queseq,
    tarseq   = tarseq,
    queries  = queries,
    offset   = offset,
    maxspace = con@alignment@dna2dna_maxspace
  )

  cds  = tarsp$CDS
  exon = tarsp$exon
  mrna = tarsp$mRNA

  hits <- gene2genome_aln$map
  met <- S4Vectors::mcols(hits)
  rng <- CNEr::second(hits)

  has_cds_match  <- GenomicRanges::findOverlaps(rng, cds)  %>% S4Vectors::queryHits() %>% unique
  has_exon_match <- GenomicRanges::findOverlaps(rng, exon) %>% S4Vectors::queryHits() %>% unique
  has_mrna_match <- GenomicRanges::findOverlaps(rng, mrna) %>% S4Vectors::queryHits() %>% unique

  met$cds_match  <- seq_along(met$score) %in% has_cds_match
  met$exon_match <- seq_along(met$score) %in% has_exon_match
  met$mrna_match <- seq_along(met$score) %in% has_mrna_match

  S4Vectors::mcols(hits) <- met

  map <-
    S4Vectors::mcols(hits) %>% as.data.frame %>%
    dplyr::mutate(target = 1:length(query)) %>%
    dplyr::select(query, target, score, logmn, pval, alnrate, cds_match, exon_match, mrna_match)

  S4Vectors::mcols(hits) <- S4Vectors::mcols(hits)[, c("strand", "query", "qwidth", "twidth")]

  gene2genome <- list(map = map, hits = hits)

  out <- list(aa2aa=aa2aa, aa2orf=aa2orf, aa2transorf=aa2transorf, gene2genome=gene2genome,
              gaps=gaps, indels=indels, summary_synder_flags=summary_synder_flags, synder_out=synder_out,
              skipped=gene2genome_aln$skipped, gene_tag=gene_tag)

  out

}


#'  Calculate significance for alignment between focal species and target genome
#'
#' @param con configuration
#' @param pair a list of alignment result
#' @export
#' @return append the list of alignment result with pvalue
calculate_match_significance <- function(pair, con){

  adjust <- function(pair, tag){
      d <- pair[[tag]]$map[, c("query", "target", "pval", "alnrate")] %>%
      dplyr::group_by(query) %>%
      dplyr::mutate(pval.adj = p.adjust(pval, method=con@alignment@padjust_method)) %>%
      dplyr::select(query, target, pval, pval.adj, alnrate)
  }

  for(tag in c("aa2aa", "aa2orf", "aa2transorf", "gene2genome")){
    pair[[sprintf("%s.pval", tag)]] <- adjust(pair=pair, tag=tag)
  }

  pair
}


#' Compare the genes from focal species to target genome and calculate significance
#'
#' @param con configuration
#' @param quesp a list of species metadata for focal species from load_species
#' @param tarsp a list of species metadata for target species from load_species
#' @param quename a string of focal species name
#' @param tarname a string of target species name
#' @export
#' @return a list of alignment result with pvalue, and save the list as rds file
secondary_data <- function(con, quesp, tarsp, quename, tarname){
  gene_tag <- load_gene_list(con@input@gene_list[[quename]])
  gene_check(quesp$mRNA, quename, gene_tag)
  pair <- compare_target_to_focal(con, quesp, tarsp, quename, tarname, gene_tag)
  pair <- calculate_match_significance(pair,con)
  saveRDS(pair,paste0(con@archive,"/",quename,"-",tarname,".rds"))
  pair
}


