#' Calls peaks downstream of GroupCoverage based on variable of interest
#'
#'
#' @param archr.proj ArchR project with group coverages
#' @param groupBy Cell population of interest which should be in cellColData
#' @param ... Additional arguments to addReproduciblePeakSet
#'
#' @return archr.proj with peaks called
#'
#' @internal
peakCalling = function(archr.proj, groupBy, ...){
 
  print(paste("Calling peaks by: ", groupBy))

  ## MACS2 peak caller
  pathToMacs2 <- ArchR::findMacs2()

  ## Define peak sets
  archr.proj = ArchR::addReproduciblePeakSet(
                ArchRProj = archr.proj, 
                groupBy = groupBy, 
                peakMethod = pathToMacs2,
                force=TRUE,
                ...)

    ## Add peak matrix
  print("Peak calling finished, adding peak matrix to ArchRProject")  
  archr.proj <- ArchR::addPeakMatrix(archr.proj)

  return(archr.proj)
}