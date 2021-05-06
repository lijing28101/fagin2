#' @import doParallel
#'
utils::globalVariables(c("%dopar%","%:%"))
NULL

#' Run fagin analysis
#'
#' @param con configuration
#' @param cores number of cpu
#' @export
#' @return a list of result for each pairwise of focal and target species and save result as excel file

run_fagin_parallel <- function(con, cores=16){

  cl <- makeCluster(cores) #not to overload your computer
  registerDoParallel(cl)

  if(!file.exists(con@archive)){
    dir.create(con@archive)
  }

  all_species <- get_species(con)

  # Step 1. load species
  foreach(species = all_species) %dopar% {

    if(!file.exists(paste0(con@archive,"/",species,"_data.rds"))){
      load_species(species,con)
    }
  }

  result <-  foreach(focal=con@input@focal_species) %:%
    foreach(target=names(con@input@syn[[focal]])) %dopar% {

      quesp <- readRDS(paste0(con@archive,"/",focal,"_data.rds"))
      tarsp <- readRDS(paste0(con@archive,"/",target,"_data.rds"))

      # Step 2. pairwise comparison
      if(!file.exists(paste0(con@archive,"/",focal,"-",target,".rds"))){

        pair <- secondary_data(con, quesp, tarsp, focal, target)

      } else {

        pair <- readRDS(paste0(con@archive,"/",focal,"-",target,".rds"))

      }

      # Step 3. summarize comparison
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
    openxlsx::saveWorkbook(wb, file = paste0(con@archive,"/",con@input@focal_species[focal],"_result.xlsx"), overwrite = TRUE)

    saveRDS(result_focal,paste0(con@archive,"/",con@input@focal_species[focal],"_result.rds"))

  }

  final_result
  stopCluster(cl)

}
