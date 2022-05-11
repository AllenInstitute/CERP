#' Creates coverages for a variable of interest  
#'
#' 
#'
#' @param ArchR_project ArchR project 
#' 
#'
#' @return ArchR_project with coverages calculated for variable of interest
#'
#' @export
GroupCoverage<-function(archr.proj = NULL, groupBy = NULL, ...){
  
  print(paste("Creating Group Coverages by: ", groupBy))
  
  archr.proj <- addGroupCoverages(
    ArchRProj = archr.proj,
    groupBy = groupBy,
    useLabels = TRUE,
    minCells = 50,
    maxCells = 10000,
    maxFragments = 25 * 10^6,
    minReplicates = 2,
    maxReplicates = 5,
    sampleRatio = 0.8,
    kmerLength = 6,
    threads = getArchRThreads(),
    returnGroups = FALSE,
    parallelParam = NULL,
    force = TRUE,
    verbose = TRUE,
    logFile = createLogFile("addGroupCoverages")
  )
  return(archr.proj)
}
