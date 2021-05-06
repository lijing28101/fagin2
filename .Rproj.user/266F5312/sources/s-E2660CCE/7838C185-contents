#' @importFrom methods isClass new as
#' @importFrom data.tree as.Node
#' @importFrom rlang .data
#' @importFrom yaml yaml.load
#' @importFrom ape read.tree
#' @importFrom L1pack l1fit
#' @importFrom fitdistrplus fitdist
#' @importFrom glue glue
#' @importFrom readr read_table
#' @importFrom rtracklayer readGFF
#' @importFrom stats quantile median sd density p.adjust
#' @importFrom CNEr first second GRangePairs
#' @importFrom magrittr is_greater_than set_names "%>%" "%T>%" "%$%" set_rownames
#' @importFrom dplyr group_by summarize mutate select filter arrange rename
#' @importFrom Biostrings subseq vcountPattern alphabetFrequency pairwiseAlignment reverseComplement pattern subject translate writeXStringSet
#' @importFrom Rsamtools scanFa getSeq indexFa FaFile
#' @importFrom GenomicRanges ranges seqnames findOverlaps GRanges strand seqinfo
#' @importFrom GenomicFeatures transcripts cds exons makeTxDbFromGRanges cdsBy extractTranscriptSeqs exonsBy
#' @importFrom IRanges reverse IRanges
#' @importFrom synder read_synmap search flag_summary as_gff
#' @importFrom S4Vectors queryHits metadata from to
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook
#' @importFrom BiocGenerics score relist
#' @importFrom utils head tail
#'
utils::globalVariables(c("%>%", ".", "%T>%"))
NULL

#' Run fagin analysis
#'
#' @param con configuration
#' @param overwrite.result a logical variable indicating whether you want to re-summarize the comparison if result already existed in dir
#' @export
#' @return a list of result for each pairwise of focal and target species and save result as excel file
run_fagin <- function(con, overwrite.result=TRUE){

  if(!file.exists(con@archive)){
    dir.create(con@archive)
  }

  all_species <- get_species(con)

  # Step 1. load species
  for(species in all_species){
    if(!file.exists(paste0(con@archive,"/",species,"_data.rds"))){

      message("load speicies ",species)

      load_species(species,con)
    } else {

      message("skip loading ", species, ", read RDS instead")

    }
  }

  final_result <- list()
  i=0

  for(focal in con@input@focal_species){

    if(overwrite.result==TRUE | !file.exists(paste0(con@archive,"/",focal,"_result.rds"))){

      message("Working on the focal species: ", focal)

      target_species <- setdiff(all_species, focal)
      quesp <- readRDS(paste0(con@archive,"/",focal,"_data.rds"))

      feature_focal <- data.frame()
      aatab_focal <- data.frame()
      orftab_focal <- data.frame()
      transtab_focal <- data.frame()
      gentab_focal <- data.frame()

      for(target in target_species){

        tarsp <- readRDS(paste0(con@archive,"/",target,"_data.rds"))

        # Step 2. pairwise comparison
        if(!file.exists(paste0(con@archive,"/",focal,"-",target,".rds"))){

          pair <- secondary_data(con, quesp, tarsp, focal, target)

        } else {

          message("skip comparison for ", focal, " and ", target, ", read RDS instead")

          pair <- readRDS(paste0(con@archive,"/",focal,"-",target,".rds"))
        }

        # Step 3. summarize comparison
        feature <- merge_feature_table(con, pair, target)
        feature_focal <- rbind(feature_focal,feature)
        hits_info <- get_hits_info(pair,con,quesp,tarsp,target)
        aatab_focal <- rbind(aatab_focal,hits_info$aatab)
        orftab_focal <- rbind(orftab_focal,hits_info$orftab)
        transtab_focal <- rbind(transtab_focal,hits_info$transtab)
        gentab_focal <- rbind(gentab_focal,hits_info$gentab)
      }

      result_focal <- list(feature_focal=feature_focal,aatab_focal=aatab_focal,orftab_focal=orftab_focal,
                           transtab_focal=transtab_focal,gentab_focal=gentab_focal)

      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "main")
      openxlsx::addWorksheet(wb, "aa2aa")
      openxlsx::addWorksheet(wb, "aa2orf")
      openxlsx::addWorksheet(wb, "aa2transorf")
      openxlsx::addWorksheet(wb, "gene2genome")
      openxlsx::writeData(wb, 1, feature_focal)
      openxlsx::writeData(wb, 2, aatab_focal)
      openxlsx::writeData(wb, 3, orftab_focal)
      openxlsx::writeData(wb, 4, transtab_focal)
      openxlsx::writeData(wb, 5, gentab_focal)
      openxlsx::saveWorkbook(wb, file = paste0(con@archive,"/",focal,"_result.xlsx"), overwrite = TRUE)

      saveRDS(result_focal,paste0(con@archive,"/",focal,"_result.rds"))

    } else {

      message("skip focal species ", focal, ", read RDS instead")

      result_focal <- readRDS(paste0(con@archive,"/",focal,"_result.rds"))

    }

    i=i+1
    final_result[[i]] <- result_focal
    names(final_result)[i] <- focal

  }

  final_result

}
