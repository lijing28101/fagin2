#' Find interval for indel
#'
#' @param si synder out
#' @param indel.threshold number for indel threshold
#' @export
#' @return GenomicRanges object
find_indels <- function(si, indel.threshold=0.25){

  "
  Find possible indels or genes of resized size

  An indel is identified based on the following two criteria:
  1. The search interval must be bounded (l_flag == 1 && r_flag == 1)
  2. The target to query ratio must be smaller than the given threshold
  "

  met <- GenomicRanges::mcols(si)

  is_bound <- met$l_flag == 1 & met$r_flag == 1
  qt_ratio <-
    GenomicRanges::width(CNEr::first(si)) /
    GenomicRanges::width(CNEr::second(si))

  met$attr[is_bound & (qt_ratio < indel.threshold)] %>% unique

}
