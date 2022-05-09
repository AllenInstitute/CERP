#' Function to calculate and export marker peaks. Exports 
#'
#' 
#'
#' @param archr.proj ArchR object with peaks called
#' @param groupBy 
#' @param output.folder If specified, will save marker heatmap pdf to output folder location, along with the peak matrix
#'
#' @return archr.proj, statistically significant peaks, peak matrix 
#' 
#' @export 
markerPeaks = function(archr.proj = NULL, groupBy = NULL, output.folder = NULL){

  ##
  dir.create(output.folder)
  
  ##
  archr.proj = addPeakMatrix(archr.proj,
                             force = TRUE)
  
  ##
  peak.matrix = getMatrixFromProject(ArchRProj = archr.proj,
                                     useMatrix = "PeakMatrix",
                                     useSeqnames = NULL,
                                     verbose = TRUE,
                                     binarize = FALSE,
                                     threads = getArchRThreads(),
                                     logFile = createLogFile("getMatrixFromProject"))
  peak.matrix@rowRanges$peak = paste(GRangesToString(peak.matrix@rowRanges,sep = c(":","-"))) 
  
  ##
  marker_features  =  getMarkerFeatures(ArchRProj = archr.proj, 
                                        useMatrix = "PeakMatrix", 
                                        groupBy = groupBy) # column results from atac integration...bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")
  
  ## Should we filter here????
  marker_peaks = as.data.table(getMarkers(marker_features))
  
  ## 
  heatmapPeaks  = plotMarkerHeatmap(seMarker = marker_features, 
                                    cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
                                    transpose = TRUE,
                                    plotLog2FC = TRUE)

  ##
  plotPDF(heatmapPeaks, name = paste0(groupBy,"_marker_peaks.pdf"), width = 16, height = 12, ArchRProj = archr.proj, addDOC = FALSE)
  
  ##
  if(!is.null(output.folder)){
    ## ? Do we need to also save the heatmap to the output foler is is the ArchRProject Plots folder sufficient?
    paste0(output.folder,"/",groupBy,"_marker_peaks_heatmap.pdf")
    pdf(file =  "/allen/programs/celltypes/workgroups/mct-t200/marcus/Z_data_analyses/MTX2060TH/MTX2060_hm_test.pdf")
    print(draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot"))
    dev.off()
    
    ## ? I don't think we want to write out the peak matrix, lets discuss
    #saveRDS(peak.matrix,file = paste0(output.folder,"/",groupBy,"_peak_matrix"))
    saveRDS(marker.list, file = paste0(output.folder,"/",groupBy,"_marker_peaks"))
  }else{}
  
  obj_list = list(archr.proj, marker.list, peak.matrix)
  names(obj_list) = c("archr.proj","marker.list","peak.matrix","heatmap")
  
  return(obj_list)
}
