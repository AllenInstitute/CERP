#' Creates coverages based on a variable of interest  
#'
#' @param archr.proj ArchR project 
#' @param groupBy Metadata field to group cells by
#' 
#' @return archr.proj with coverages calculated for variable of interest
#'
#' @keywords internal
groupCoverage<-function(archr.proj, groupBy, ...){
  
  print(paste("Creating Group Coverages by: ", groupBy))
  
  ## Add groups to Arrow files
  archr.proj <- addGroupCoverages(
    ArchRProj = archr.proj,
    groupBy = groupBy,
    force = TRUE,
    ...
  )
  return(archr.proj)
}
