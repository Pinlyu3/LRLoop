#' Write the LRloop_info into .csv files
#' @name writeLRL
#' @param LRloop_info: The list of the LRloop_network info calculated by function "LRL_info_collection"
#' @param filedir: Directory to save the files
#' @import nichenetr tidyverse Seurat
#' @export


writeLRL <- function(LRloop_info, filedir) {
  for (i in 1:length(LRloop_info)) {
    chart = LRloop_info[[i]]
    name = names(LRloop_info)[i]
    write.csv(chart, file = sprintf("%s/%s.csv", filedir, name), row.names = FALSE)
  }
}