#' Function to calculate gini index - returns marker peaks file with gini_index as a new column
#' (Marcus Hooper)
#' @param archr.proj ArchR object with peaks called
#' @param groupBy Metadata field to group cells by
#' @param archr.visualize Should marker peak visualizations be produced
#' @param filename # output file for markerpeaks
#' @return marker peak table


library(ArchR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(compiler)
library(parallel)
library(data.table)
library(DescTools)
library(doParallel)

#main function for returning markerpeaks with gini indexes
calculate_gini_index = function(archr.proj, 
                                groupBy,
                                filename = file.path(file.path(getOutputDirectory(archr.proj), "MarkerPeaks", groupBy), paste0(groupBy, "_annotated_markerPeaks.tsv")),
                                max_sample_size = 500){

  # read marker peaks -------------------------------------------------------
  marker_peaks = read.table(filename,sep = "\t",header = T)
  # subsample cells ---------------------------------------------------------
  print("sampling cells")
  data = getCellColData(archr.proj)
  data$cellNames = rownames(data)
  data = as.data.table(data)
  data[,`:=` (cells_per_group = .N), by = groupBy]
  data$groupBy = data[,groupBy, with = F]
  cell_count_data = unique(data[,c('cells_per_group', 'groupBy'), with = F])
  #get # of samples per group, use group size if smaller than max_sample_size
  cell_count_data$sample_size = ifelse(max_sample_size < cell_count_data$cells_per_group, max_sample_size,no =  cell_count_data$cells_per_group)
  cellNames_for_calculation = unlist(lapply(FUN = sample_cells,X = unique(data$groupBy),cell_count_data = cell_count_data, data = data))
  arch = ArchR::subsetCells(ArchRProj = archr.proj, cellNames = cellNames_for_calculation)
  
  # get peak matrix ---------------------------------------------------------
  print("getting peak matrix")
  peak_matrix <-getMatrixFromProject(ArchRProj = arch,useMatrix = "PeakMatrix",useSeqnames = NULL,verbose = TRUE,binarize = FALSE,threads = getArchRThreads(),logFile = createLogFile("getMatrixFromProject"))
  peak_matrix@rowRanges$peak = paste0(as.character(peak_matrix@rowRanges@seqnames),":",as.character(peak_matrix@rowRanges@ranges))
  peak_matrix = peak_matrix[peak_matrix@rowRanges$peak %in% unique(marker_peaks$coordinates)]
  peak_matrix@colData$cellNames<-colnames(peak_matrix)
  print("peak_matrix_generated")
  
  n_cells = length(colnames(peak_matrix))
  n_groups = length(unique(as.character(peak_matrix@colData[groupBy][,1])))
  

  #get gini indexes parallel ----------------------------------------------
  numWorkers <- 10 # may change this later, will run out of memory if this is too high
  cl <-makeCluster(numWorkers, type="PSOCK")
  clusterCall(cl, function() library(DescTools))
  clusterExport(cl=cl, varlist=c("peak_matrix", "n_cells", "n_groups","groupBy"),envir=environment())
  clusterEvalQ(cl, library("data.table"))
  print("calculating gini index")
  time = Sys.time()
  
  peak_indices = 1:length(peak_matrix@rowRanges$peak)
  gini_indexes = parLapply(cl = cl,X = peak_indices,fun = get_gini_index_peak, peak_matrix = peak_matrix)

  time_end = Sys.time()
  print('time to calculate gini index:')
  print(time_end - time)
  stopCluster(cl)
  
  gini_indexes = data.table(rbindlist(gini_indexes))
  marker_peaks = as.data.table(marker_peaks)
  gini_indexes$coordinates = gini_indexes$peak
  marker_peaks = merge.data.table(marker_peaks, gini_indexes[,.(coordinates,peak_signal_total,gini_index)], by = "coordinates")
  
  if(!is.null(filename)){
    write.table(marker_peaks, row.names=F, col.names=T, sep="\t", file=filename)
  }
}

#subsamples cells for gini index calculation
sample_cells = function(group,cell_count_data =NULL, data =NULL){
  sample_size = cell_count_data[groupBy == group,]$sample_size
  return(sample(data[groupBy == group,]$cellNames,sample_size))
}

#function gets gini index for given peak
get_gini_index_peak = function(peak_index,peak_matrix){
  archR_sparse<-peak_matrix ##loaded from above ##peak data is the result of saving archR data (getting peak matrix)
  peak_signal = data.table(data.frame("peak" = rowRanges(archR_sparse[peak_index,])$peak,"peak_signal" = matrix(archR_sparse[peak_index,]@assays@data@listData$PeakMatrix),"groupBy" = as.character(archR_sparse@colData[groupBy][,1]),row.names = colnames(archR_sparse)),keep.rownames = T)
  colnames(peak_signal) = c("cellNames","peak","peak_signal","cell.population")
  peak_signal[,`:=` ( peak_signal_total = sum(peak_signal) ), by = 'peak']
  peak_signal[, `:=` ( peak_signal_group_total = sum(peak_signal)),by = c("cell.population",'peak')]
  peak_signal[, `:=` ( peak_signal_group = mean(peak_signal)),by = c("cell.population",'peak')]
  peak_signal[, `:=` ( cells_per_group = .N),by = c("cell.population")]
  peak_signal[, `:=` ( cells_total = n_cells)]
  peak_signal$cell_prop = peak_signal$cells_per_group / peak_signal$cells_total
  peak_signal_by_group = unique(peak_signal[,.(peak,peak_signal_total,peak_signal_group,cell.population,cell_prop,peak,cells_per_group)])
  
  peak_signal_by_group$gini_index = rep(Gini(x = peak_signal_by_group$peak_signal_group, n = peak_signal_by_group$cells_per_group),length(peak_signal_by_group$cell.population))
  return(peak_signal_by_group[1,.(peak,peak_signal_total,gini_index)])
}
get_gini_index_peak = cmpfun(get_gini_index_peak)
