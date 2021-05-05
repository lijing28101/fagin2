#' Make SequenceInfo
#'
#' @param x XStringSet
#' @param species_name string
#' @return SequenceInfo
make_seqinfo <- function(x, species_name){
  si <- GenomeInfoDb::seqinfo(x)
  GenomeInfoDb::genome(si) <- species_name
  si
}

#' Filter uname entries
#'
#' @param gr GenomicRanges object
#' @param label string
#' @return GenomicRanges object
filter_with_warning__unnamed_entries <- function(gr, label=NULL){

  "Check for unnamed entries in a GRanges object. If there are any, remove them
  and raise a warning."

  if(any(is.na(names(gr)))){
    n_cds_with_unnamed_mRNA <- gr[is.na(names(gr))] %>% length
    total <- GenomicRanges::seqnames(gr) %>% length
    warning(glue::glue(.sep=" ",
      "{.label(label)}: {n_cds_with_unnamed_mRNA} out of {total} entries",
      "have no name associated with them. This may be bad.",
      "All of these entries will be removed from the analysis"
    ))
    gr <- gr[!is.na(names(gr))]
  }

  gr
}


#' Convert fasta file to XStringSet object
#'
#' @param x indexed FASTA file
#' @return XStringSet object
convert_FaFile_to_XStringSet <- function(x, ...){

  "FaFile -> XStringSet

  Load an XStringSet object from an indexed FASTA file. Use the first word in
  the header as the sequence name

  A wrapper for scanFa that sets the seqname to the first word in the header
  This is probably the behaviour the function should have, since this is done
  when subsetting a sequence using a GRanges object.
  "

  seq <- Rsamtools::scanFa(x, ...)
  names(seq) <- sub(" .*", "", names(seq))
  seq
}

#' Get transcripts sequence from GenomicRanges object
#'
#' @param x indexed FaFile
#' @param species_name string
#' @param transcripts GenomicRanges object
#' @return indexed FaFile
get_trans_dna <- function(x, species_name, transcripts){
  GenomicFeatures::extractTranscriptSeqs(x, transcripts) %>%
    {

      "Ensure all transcripts are named. If any of the names are missing (NA),
    then the Biostrings::writeXStringSet will report an error as a note,
    which won't stop processing. By checking here I can stop analysis at the
    right time.

    If any transcripts do have missing names, then these are removed with a
    warning.
    "

      na_indices <- which(is.na(names(.)))

      if(length(na_indices) > 0){ msg <-
        "%s of %s transcripts had missing names. This is not good. You should look into
the problem. For now, the offending transcripts have been removed."
      warning(sprintf(msg, length(na_indices), length(.)))
      . <- .[-na_indices]
      }

      .

    } %>% {

      "Print the transcripts to a temporary file"

      filepath <- file.path(con@archive, paste0(".", species_name, "_trans.fna"))
      Biostrings::writeXStringSet(., filepath=filepath)
      filepath

    } %>% {

      "Build an indexed version of the genome"

      Rsamtools::indexFa(.)

      .

    } %>% {

      "Get a reference to the genome"

      Rsamtools::FaFile(.)

    }
}
