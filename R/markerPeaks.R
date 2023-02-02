#' Function to calculate and export marker peaks. Exports 
#'
#' @param archr.proj ArchR object with peaks called
#' @param groupBy Metadata field to group cells by
#' @param unique.id Unique peak calling run ID
#'
#' @return marker peak table
#' 
#' @keywords internal
markerPeaks = function(archr.proj, groupBy, unique.id=NULL){

  ## Stat. test for marker peaks
  marker_features  =  getMarkerFeatures(ArchRProj = archr.proj, 
                                        useMatrix = "PeakMatrix", 
                                        groupBy = groupBy)

  ## Extract marker peaks 
  marker.list = as.data.table(getMarkers(marker_features))

  ## add genome_coordinates
  marker.list$genome_coordinates = paste0(marker.list$seqnames, ":" , marker.list$start, "-" , marker.list$end)

  ## Save marker peaks
  write.table(marker.list, 
              file=file.path(getOutputDirectory(archr.proj), "MarkerPeaks", groupBy, paste0(groupBy, "_markerPeaks.tsv")), 
              sep="\t", 
              row.names=FALSE)

  ## Record marker list
  archr.proj@projectMetadata[[paste0(groupBy, "_markerPeaks")]] = marker.list

  return(archr.proj)
}
