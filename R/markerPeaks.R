#' Function to calculate and export marker peaks. Exports 
#'
#' @param archr.proj ArchR object with peaks called
#' @param groupBy Metadata field to group cells by
#' @param archr.visualize Should marker peak visualizations be produced
#' @param output.dir ArchR directory in which to save marker peak file
#'
#' @return archr.proj, statistically significant peaks, peak matrix 
#' 
#' @keywords internal
markerPeaks = function(archr.proj, groupBy, archr.visualize=TRUE, output.dir="MarkerPeaks"){

  ## Stat. test for marker peaks
  marker_features  =  getMarkerFeatures(ArchRProj = archr.proj, 
                                        useMatrix = "PeakMatrix", 
                                        groupBy = groupBy)

  ## Extract marker peaks 
  marker.list = as.data.table(getMarkers(marker_features))
  
  ## add genome_coordinates
  marker.list$genome_coordinates = paste0(marker.list$seqnames,":",marker.list$start,"-",marker.list$end)
  
  ## Plot marker peaks
  if(archr.visualize){
    heatmapPeaks  = plotMarkerHeatmap(seMarker = marker_features, 
                                      cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
                                      transpose = TRUE,
                                      plotLog2FC = TRUE)

    ##
    plotPDF(heatmapPeaks, name = paste0(groupBy,"_marker_peaks.pdf"), width = 16, height = 12, ArchRProj = archr.proj, addDOC = FALSE)
  }

  ## Create new directory under ArchR project for marker peak results
  dir.create(file.path(getOutputDirectory(archr.proj), output.dir), showWarnings = FALSE)
  dir.create(file.path(getOutputDirectory(archr.proj), output.dir, groupBy), showWarnings = FALSE)

  ## Save marker peaks
  marker.file = file.path(getOutputDirectory(archr.proj), output.dir, groupBy, paste0(groupBy, "_markerPeaks.tsv"))
  write.table(marker.list, file=marker.file, sep="\t", row.names=FALSE)
}
