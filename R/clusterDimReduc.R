#' Performs LSI, clustering and UMAP
#'
#' Sets global variables about specific genome build and number of threads to utilze on the machine.
#'
#' @param archr.proj ArchRProject object
#' @param varFeatures Number of variable features to use for cluster
#' @param resolution Vector of clustering resolutions for multi-round clusters
#' @param ... Additional arguments for addIterativeLSI
#' 
#' @return None
#'
#' @keywords internal
clusteringDimReduc = function(archr.proj, varFeatures=15000, resolution=c(0.2, 1, 2), ...){

  ## Perform clustering
  archr.proj = addIterativeLSI(
    ArchRProj = archr.proj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = length(resolution), 
    clusterParams = list( 
      resolution = resolution, 
      sampleCells = 10000, 
      n.start = 10
    ), 
    varFeatures = varFeatures, 
    force = TRUE,
    ...
  )

  ## Cluster
  archr.proj = addClusters(input = archr.proj, reducedDims = "IterativeLSI", force=TRUE)

  ## Compute low dimensional embedding
  archr.proj = addUMAP(ArchRProj = archr.proj, reducedDims = "IterativeLSI", force=TRUE)

  return(archr.proj)
}