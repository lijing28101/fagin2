#' Get alignment info
#'
#' @param xs PairwiseAlignments
#' @export
#' @return data frame for alignment info
get_alignment_table <- function(xs){
  data.frame(
    seqid = xs@metadata$query,
    target = xs@metadata$target,
    qstart = xs@pattern@range@start,
    qwidth = xs@pattern@range@width,
    qseq = as.character(xs@pattern),
    tstart = xs@subject@range@start,
    twidth = xs@subject@range@width,
    tseq = as.character(xs@subject),
    score = xs@score
  )
}


#' Get query gene position
#'
#' @param dat synder GFF object
#' @export
#' @return data frame for query gene position
get_query_position <- function(dat){
  a <- as.data.frame(dat)
  data.frame(
    seqid=a$attr,
    qchrom=a$seqnames,
    qlocstart=a$start,
    qlocstop=a$end,
    qstrand=a$strand
  )
}


#' Get target gene/ORF position
#'
#' @param dat synder GFF object
#' @export
#' @return data frame for target gene/ORF position
get_target_position <- function(dat){
  a <- as.data.frame(dat)
  data.frame(
    target=a$attr,
    tchrom=a$seqnames,
    tlocstart=a$start,
    tlocstop=a$end,
    tstrand=a$strand
  )
}


#' Get target mRNA position
#'
#' @param dat synder GFF object
#' @export
#' @return data frame for target mRNA position
get_target_mRNA_position <- function(dat){
  a <- as.data.frame(dat)
  data.frame(
    target=a$attr,
    target_trans=a$seqnames,
    tlocstart=a$start,
    tlocstop=a$end,
    tstrand=a$strand
  )
}


#' Summary for the gene with significant hits
#'
#' @param pair a list of alignment result
#' @param con configuration
#' @param quesp a list of species metadata for focal species from load_species
#' @param tarsp a list of species metadata for target species from load_species
#' @param tarname a string of target species name
#' @export
#' @return a list of summary for the gene with significant hits
get_hits_info <- function(pair,con,quesp,tarsp,tarname){

  p2p_cutoff <- con@alignment@pvalue@prot2prot
  p2a_cutoff <- con@alignment@pvalue@prot2allorf
  d2d_cutoff <- con@alignment@pvalue@dna2dna
  p2t_cutoff <- con@alignment@pvalue@prot2transorf

  p2p_alncut <- con@alignment@alnrate@prot2prot
  p2a_alncut <- con@alignment@alnrate@prot2allorf
  d2d_alncut <- con@alignment@alnrate@dna2dna
  p2t_alncut <- con@alignment@alnrate@prot2transorf

  aatab <- pair[["aa2aa.pval"]] %>%
    dplyr::filter(pval.adj < p2p_cutoff & alnrate > p2p_alncut) %>%
    dplyr::select(seqid=query, target=target, pval.adj, alnrate) %>%
    merge(., get_alignment_table(pair$aa2aa$aln), by=c("seqid", "target")) %>%
    dplyr::mutate(target_species=tarname)

  orftab <- pair[["aa2orf.pval"]] %>%
    dplyr::filter(pval.adj < p2a_cutoff & alnrate > p2a_alncut) %>%
    dplyr::select(seqid=query, target=target, pval.adj, alnrate) %>%
    merge(., get_alignment_table(pair$aa2orf$aln), by=c("seqid", "target")) %>%
    merge(., get_query_position(quesp[["mRNA"]]), by="seqid", all.x=T) %>%
    merge(., get_target_position(tarsp[["orfgff"]]), by="target", all.x=T) %>%
    dplyr::select(seqid,qchrom,qlocstart,qlocstop,qstrand,qstart,qwidth,target,tchrom,tlocstart,tlocstop,tstrand,tstart,twidth,score,pval.adj,alnrate,qseq,tseq) %>%
    unique() %>%
    dplyr::mutate(target_species=tarname)

  target_transoff_loc <- get_target_mRNA_position(tarsp[["transorfgff"]]) %>%
    merge(.,dplyr::rename(get_query_position(tarsp[["mRNA"]]),target_trans=seqid),by="target_trans",all.x=T) %>%
    select(target_trans,target,tchrom=qchrom,tlocstart,tlocstop,tstrand=qstrand)

  transtab <- pair[["aa2transorf.pval"]] %>%
    dplyr::filter(pval.adj < p2t_cutoff & alnrate > p2t_alncut) %>%
    dplyr::select(seqid=query, target=target, pval.adj, alnrate) %>%
    merge(., get_alignment_table(pair$aa2transorf$aln), by=c("seqid", "target")) %>%
    merge(., get_query_position(quesp[["mRNA"]]), by="seqid", all.x=T) %>%
    merge(., target_transoff_loc, by="target", all.x=T) %>%
    dplyr::select(seqid,qchrom,qlocstart,qlocstop,qstrand,qstart,qwidth,target,target_trans,tchrom,tlocstart,tlocstop,tstrand,tstart,twidth,score,pval.adj,alnrate,qseq,tseq) %>%
    unique() %>%
    dplyr::mutate(target_species=tarname)

  map <- pair[["gene2genome"]][["map"]]
  gentab <- as.data.frame(pair[["gene2genome"]][["hits"]][map$target]) %>%
    dplyr::select(
      seqid  = query,
      target = second.seqnames,
      tstart = second.start,
      twidth = second.width,
      qstart = first.start,
      qwidth = first.width,
      strand = strand
    ) %>%
    dplyr::mutate(id=map$target) %>%
    merge(.,pair[["gene2genome.pval"]] %>%
            dplyr::select(seqid=query, id=target, pval.adj, alnrate),
          by=c("seqid", "id")
    ) %>%
    dplyr::filter(pval.adj < d2d_cutoff & alnrate > d2d_alncut) %>%
    dplyr::select(seqid,qstart,qwidth,strand,target,id,tstart,twidth,pval.adj,alnrate) %>%
    unique() %>%
    dplyr::mutate(target_species=tarname)

  hits_info <- list(aatab=aatab,orftab=orftab,transtab=transtab,gentab=gentab)
  hits_info

}
