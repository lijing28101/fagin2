
#' Do alignment for new genes from focal species agaginst to target genome
#'
#' @param con configuration
#' @param quesp a list of species metadata for focal species from load_species
#' @param tarsp a list of species metadata for target species from load_species
#' @param quename a string of focal species name
#' @param tarname a string of target species name
#' @param gene_tag a list of genes as character vector
#' @param pair existing pair
#' @export
#' @return a list of alignment result
add_gene <- function(con, quesp, tarsp, quename, tarname, gene_tag, pair){

  message("Processing ", quename, " vs ", tarname)

  nsims_prot <- con@alignment@simulation@prot2prot
  nsims_allorf <- con@alignment@simulation@prot2allorf
  nsims_transorf <- con@alignment@simulation@prot2transorf
  subMat <- con@alignment@substitution_matrix

  # SyntenyMap -> mRNA -> SearchIntervals
  synder_out <- pair$synder_out


  f_si_map <- pair$f_si_map

  f_si_map_orf <- pair$f_si_map_orf

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

  f_si_map_transorf <- pair$f_si_map_transorf

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
  tarseq <- pair$tarseq
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


  pair$aa2aa=aa2aa
  pair$aa2orf=aa2orf
  pair$aa2transorf=aa2transorf
  pair$gene2genome=gene2genome
  pair$skipped=gene2genome_aln$skipped
  pair$gene_tag=gene_tag

  pair
}


#' Compare the new genes from focal species to target genome and calculate significance
#'
#' @param con configuration
#' @param quesp a list of species metadata for focal species from load_species
#' @param tarsp a list of species metadata for target species from load_species
#' @param quename a string of focal species name
#' @param tarname a string of target species name
#' @param pair existing pair
#' @param gene_path path to the gene list
#' @export
#' @return a list of alignment result with pvalue, and save the list as rds file
secondary_data_add <- function(con, quesp, tarsp, quename, tarname, pair, gene_path){
  gene_tag <- load_gene_list(gene_path)
  gene_check(quesp$mRNA, quename, gene_tag)
  pair <- add_gene(con, quesp, tarsp, quename, tarname, gene_tag, pair)
  pair <- calculate_match_significance(pair,con)
  pair$feature <- merge_feature_table(con, pair, target)

  #saveRDS(pair,paste0(con@archive,"/",prefix,"/",quename,"-",tarname,".rds"))
  pair
}


#' Add new gene list to fagin analysis
#'
#' @param con configuration
#' @param prefix path to dir for saving new result
#' @param gene_path path to the gene list
#' @param cores number of cpu to be used
#' @param cl.type cluster type, default is "FORK", Windows should set "PSOCK"
#' @export
#' @return a list of result for each pairwise of focal and target species and save result as excel file
run_fagin_add <- function(con, gene_path, prefix, cores=16,cl.type="FORK"){

  cl <- makeCluster(cores, type = cl.type) #not to overload your computer
  registerDoParallel(cl)


  all_species <- get_species(con)

  result <-  foreach(focal=con@input@focal_species) %:%
    foreach(target=names(con@input@syn[[focal]])) %dopar% {

      .GlobalEnv$con = con


        quesp <- readRDS(paste0(con@archive,"/",focal,"_data.rds"))
        tarsp <- readRDS(paste0(con@archive,"/",target,"_data.rds"))
        pair <- readRDS(paste0(con@archive,"/",focal,"-",target,".rds"))

        pair <- secondary_data_add(con, quesp, tarsp, focal, target, pair, gene_path)


      feature <- merge_feature_table(con, pair, target)

      hits_info <- get_hits_info(pair,con,quesp,tarsp,target)
      hits_info[["feature"]] <- feature
      hits_info
    }

  listname <- c("aatab","orftab","transtab","gentab","feature")

  final_result <- list()

  for(focal in 1:length(con@input@focal_species)) {
    result_focal <- list()
    for(i in listname){
      tab <- lapply(result[[focal]], function(x) x[[i]])
      tabdf <- do.call("rbind", tab)
      result_focal[[i]] <- tabdf
    }
    final_result[[con@input@focal_species[focal]]] <- result_focal

    if(!file.exists(paste0(con@archive,"/",prefix))){
      dir.create(paste0(con@archive,"/",prefix))
    }

    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "main")
    openxlsx::addWorksheet(wb, "aa2aa")
    openxlsx::addWorksheet(wb, "aa2orf")
    openxlsx::addWorksheet(wb, "aa2transorf")
    openxlsx::addWorksheet(wb, "gene2genome")
    openxlsx::writeData(wb, 1, result_focal[["feature"]])
    openxlsx::writeData(wb, 2, result_focal[["aatab"]])
    openxlsx::writeData(wb, 3, result_focal[["orftab"]])
    openxlsx::writeData(wb, 4, result_focal[["transtab"]])
    openxlsx::writeData(wb, 5, result_focal[["gentab"]])
    openxlsx::saveWorkbook(wb, file = paste0(con@archive,"/",prefix,"/",con@input@focal_species[focal],"_result.xlsx"), overwrite = TRUE)

    saveRDS(result_focal,paste0(con@archive,"/",prefix,"/",con@input@focal_species[focal],"_result.rds"))
  }

  stopCluster(cl)
  return(final_result)

}


