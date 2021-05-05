#' Find all features in the gff GenomicRanges object that overlap a search interval
#'
#' @param gff GenomicRanges object
#' @param si synder out
#' @param type a string from any, start, end, within, equal
#' @export
#' @return data frame of overlap features
overlapMap <- function(gff, si, type=NULL){

  "Find all features in the gff GenomicRanges object that overlap a search interval."

  if(! class(gff) == "GFF"){
    stop("Expected 'GFF', got '", class(gff), "'")
  }

  # A Hits object
  # from(o) accesses the si indices
  # to(o) acceses the ft indices
  GenomicRanges::findOverlaps(
      query = CNEr::second(si)
    , subject = gff
    , type = "any" # from: [any, start, end, within, equal]
    , select = "all" # from: [all, first, last, arbitrary]
    , ignore.strand = TRUE # This MUST be set to TRUE. The search intervals are
                           # always '+', which is arbitrary since strand does
                           # not matter for them. If strand is considered here,
                           # all minus strand mRNAs will be ignored.
  ) %>% {
    data.frame(
      query            = GenomicRanges::mcols(si[S4Vectors::from(.)])$attr,
      target           = GenomicRanges::mcols(gff[S4Vectors::to(.)])$attr,
      type             = GenomicRanges::mcols(gff[S4Vectors::to(.)])$type,
      qid              = S4Vectors::from(.),
      tid              = S4Vectors::to(.),
      stringsAsFactors = FALSE
    )
  }
}
