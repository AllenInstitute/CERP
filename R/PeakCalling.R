#' Calls peaks downstream of GroupCoverage based on variable of interest
#'
#' 
#'
#' @param ArchR_project ArchR project with group coverages
#' @param cell_population cell population of interest this should be in cellColData
#'
#' @return ArchR_project with peaks called
#'
#' @export
library(ArchR)

PeakCalling<-function(ArchR_project = NULL,cell_population = NULL){
  dir.create(output_folder)
 
  print(paste("Calling peaks by",cell_population))
  ArchR_project <- addReproduciblePeakSet(ArchRProj = ArchR_project, groupBy = cell_population, peakMethod = "Macs2",minCells = 50,force = T)
  print("Peak calling finished")
  
  
  return(ArchR_project)
  
}