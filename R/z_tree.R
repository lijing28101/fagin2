#' Get all species in the study
#'
#' @param con configuration
#' @export
#' @return vector of all species
get_species <- function(con){
  ape::read.tree(con@input@tree)$tip.label %>%
    gsub(pattern=" ", replacement="_")
}
