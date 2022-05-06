#' Function to calculate and export marker peaks. Exports 
#'
#' 
#'
#' @param ArchR_project ArchR object with peaks called
#' @param output_folder If specified, will save marker heatmap pdf to output folder location, along with the peak matrix
#'
#' @return list, ArchR_project, statistically significant peaks, peak matrix 
#' 
#' @export 
library(ArchR)
library(data.table)
library(Signac)


preprocess_by_var<-function(ArchR_project = NULL,cell_population = NULL,output_folder = NULL){
  dir.create(output_folder)
  

  ArchR_project<-addPeakMatrix(ArchR_project,force = T)
  
  peak_matrix<-getMatrixFromProject(ArchRProj = ArchR_project,useMatrix = "PeakMatrix",
                                    useSeqnames = NULL,verbose = TRUE,binarize = FALSE,threads = getArchRThreads(),
                                    logFile = createLogFile("getMatrixFromProject"))
  peak_matrix@rowRanges$peak<-paste(GRangesToString(peak_matrix@rowRanges,sep = c(":","-"))) 
  
  #saveRDS(peak_matrix, file = paste0(output_folder,'/',cell_population,'_peakmatrix'))
  marker_features <- getMarkerFeatures(ArchRProj = ArchR_project, useMatrix = "PeakMatrix", groupBy = cell_population) # column results from atac integration...bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")
  marker_peaks<-as.data.table(getMarkers(marker_features))
  

  heatmapPeaks <-plotMarkerHeatmap(seMarker = marker_features, cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
                                    transpose = TRUE,plotLog2FC = T)
  
  draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  
  plotPDF(heatmapPeaks, name = paste0(cell_population,"_marker_peaks.pdf"), width = 16, height = 12, ArchRProj = ArchR_project, addDOC = FALSE)
  hm<-draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

  if(!is.null(output_folder)){
    paste0(output_folder,"/",cell_population,"_marker_peaks_heatmap.pdf")
    pdf(file =  "/allen/programs/celltypes/workgroups/mct-t200/marcus/Z_data_analyses/MTX2060TH/MTX2060_hm_test.pdf")
    hm
    dev.off()
    saveRDS(peak_matrix,file = paste0(output_folder,"/",cell_population,"_peak_matrix"))
    saveRDS(marker_list,file = paste0(output_folder,"/",cell_population,"_marker_peaks"))
  }else{}
  
  obj_list<-list(ArchR_project, marker_list,peak_matrix,hm)
  names(obj_list)<-c("ArchR_project","marker_list","peak_matrix","heatmap")
  
  return(obj_list)
  
}
