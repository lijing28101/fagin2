#' Summarize functions
#'
#' @param x An object to summarized
#' @name fagin_summary
NULL


#' @rdname fagin_summary
#' @export
summarize_numeric <- function(x){
  if(length(x) == 0){
    return(NULL)
  }
  new(
    "numeric_summary",
    min              = min(x),
    q25              = stats::quantile(x, probs=0.25),
    median           = stats::median(x),
    q75              = stats::quantile(x, probs=0.75),
    max              = max(x),
    mean             = mean(x),
    sd               = stats::sd(x),
    n                = length(x),
    density          = stats::density(x, kernel="gaussian")
  )
}

#' @rdname fagin_summary
#' @export
summarize_faa <- function(x){
  if(length(x) == 0){
    return(faa_summary())
  }

  table <- data.frame(
    seqids = names(x),
    length = Biostrings::width(x)
  )

  new(
    "faa_summary",
    initial_residue = Biostrings::subseq(x, start=1,  width=1) %>%
                      as.character %>% as.factor %>% summary(maxsum=Inf),
    final_residue = Biostrings::subseq(x, start=-1, width=1) %>%
                    as.character %>% as.factor %>% summary(maxsum=Inf),
    has_internal_stop = Biostrings::subseq(x, start=1, end=-2) %>%
                        Biostrings::vcountPattern(pattern='*') %>%
                        magrittr::is_greater_than(0),
    # inherited from seq_summary
    table = table,
    comp = Biostrings::alphabetFrequency(x)
  )
}



#' @rdname fagin_summary
#' @export
summarize_dna <- function(x){

  if(class(x) == 'FaFile')
    x <- Rsamtools::scanFa(x)

  if(length(x) == 0)
    return(dna_summary())

  table <- data.frame(
    seqids = names(x),
    length = Biostrings::width(x)
  )

  new(
    "dna_summary",
    n_triple      = (Biostrings::width(x) %% 3 == 0) %>% sum,
    initial_codon = Biostrings::subseq(x, start=1,  width=3) %>%
                    as.character %>% as.factor %>% summary(maxsum=Inf),
    final_codon   = Biostrings::subseq(x, start=-3,  width=3) %>%
                    as.character %>% as.factor %>% summary(maxsum=Inf),
    # inherited from seq_summary
    table         = table,
    comp          = Biostrings::alphabetFrequency(x)
  )
}


#' @rdname fagin_summary
#' @export
summarize_phase <- function(phases, aa){
  # This should always be true. If it is not, there is a logical problem
  # in the code, not the data.
  stopifnot(length(aa) == length(phases))
  new(
    "phase_summary",
    table = factor(phases) %>% table,
    incomplete_models = names(aa)[phases != 0]
  )
}

#' @rdname fagin_summary
#' @export
summarize_granges <- function(x){

  if(length(x) == 0){
    return(granges_summary())
  }

  xdf <- data.frame(
    seqid = x@seqnames,
    stop  = GenomicRanges::end(x),
    start = GenomicRanges::start(x)
  )

  table <- dplyr::group_by(xdf, .data$seqid) %>% 
    dplyr::summarize(
      min   = min(.data$start),
      max   = max(.data$stop)
    )

  new(
    "granges_summary",
    table = table,
    width = summarize_numeric(GenomicRanges::width(x))
  )
}

#' @rdname fagin_summary
#' @export
summarize_gff <- function(x){

  feat_trans <- GenomicFeatures::transcripts(x)
  feat_cds <- GenomicFeatures::cds(x)
  feat_exons <- GenomicFeatures::exons(x)

  table <-
    GenomicRanges::as.data.frame(feat_trans) %>%
    dplyr::group_by(.data$seqnames) %>%
    dplyr::summarize(
      min   = min(.data$start),
      max   = max(.data$end),
      mRNAs = length(.data$start)
    )

  new(
    "gff_summary",
    table       = table,
    mRNA_length = summarize_numeric(feat_trans %>% IRanges::width()),
    CDS_length  = summarize_numeric(feat_exons %>% IRanges::width()),
    exon_length = summarize_numeric(feat_cds   %>% IRanges::width())
  )
}

#' @rdname fagin_summary
#' @export
summarize_nstring <- function(x){
  # TODO: should pass more info, need a new numeric summary class
  if(class(x) == "GRanges"){
    x <- GenomicRanges::ranges(x)
  }
  if(length(x) > 0) {
    s <- summary(IRanges::width(x))
  } else {
    s <- summary(NA_real_)
  }
  ss <- as.vector(s)
  names(ss) <- names(s)
  ss
}

#' @rdname fagin_summary
#' @export
summarize_syn <- function(x){

  stopifnot(class(x) == "Synmap")

  if(length(x) == 0)
    return(synmap_summary())

  qwidth <- GenomicRanges::width(CNEr::first(x))
  twidth <- GenomicRanges::width(CNEr::second(x))
  score  <- GenomicRanges::mcols(x)$score
  new(
    "synmap_summary",
    nrow                    = length(x),
    width                   = summarize_numeric(qwidth),
    score                   = summarize_numeric(score),
    query_target_log2_ratio = summarize_numeric(log2 (qwidth / twidth) )
  )
}
