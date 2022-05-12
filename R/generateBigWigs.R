#' Performs LSI, clustering and UMAP
#'
#' Sets global variables about specific genome build and number of threads to utilze on the machine.
#'
#' @param archr.proj ArchRProject object
#' @param groupBy Metadata field to group cells by
#' @param tileSize Width of tiles that will section the genome
#' @param normMethod Normalization method for fragment pileups 
#' @param maxCells Number of cells to use for each group.
#' 
#' @return None
#'
#' @keywords internal
clusteringDimReduc = function(archr.proj, groupBy, tileSize=25, normMethod="ReadsInTSS", maxCells=NULL){

  ## Build BigWig files post group coverage
  getGroupBW(
    ArchRProj = archr.proj,
    groupBy = groupBy,
    normMethod = normMethod,
    tileSize = tileSize,
    maxCells = maxCells,
    ceiling = 4,
    verbose = TRUE,
    threads = getArchRThreads(),
    logFile = createLogFile("getGroupBW")
  )

  ## Generate fragment files

}