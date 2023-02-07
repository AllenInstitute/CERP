#'gets peak to gene links and merge with existing marker peaks'
# archr.proj must have ar_id and ar_directory for this to work properly. Future updates will allow for counts matrices to generate gene expression matrices.
#' @param archr.proj 
#' @param max_distance distance from the gene to search for linked peaks - default 250kb
#' @param reduced_dims reduced dimension parameter for p2gene links ArchR function, default "LSI_Combined"
#' @param useMatrix' matrix to use for constructing p2gene links, can be GeneExpressionMatrix (default) or GeneIntegrationMatrix
#' @param marker.file' annotated marker peaks file path

peak2Gene = function(archr.proj = NULL, max_distance = 250000, reduced_dims = "LSI_Combined", useMatrix = "GeneExpressionMatrix",
                     marker.file = file.path(file.path(getOutputDirectory(archr.proj), "MarkerPeaks", groupBy), paste0(groupBy, "_annotated_markerPeaks.tsv"))){
  #library(magrittr)
  #load markerPeaks
  #check for LSI if not there, run LSI 
  archr.proj <- addPeak2GeneLinks(
    ArchRProj = archr.proj,
    reducedDims = reduced_dims,
    useMatrix = useMatrix,
    maxDist = max_distance,
  )
  p2geneDF <- metadata(archr.proj@peakSet)$Peak2GeneLinks
  p2geneDF$correlated_gene <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
  p2geneDF$coordinates <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), ":", start(.), "-", end(.))})[p2geneDF$idxATAC]
  p2geneDF$cor_FDR<-p2geneDF$FDR
  p2geneDF = as.data.table(p2geneDF)
  p2geneDF = p2geneDF[,.(coordinates,correlated_gene,Correlation,cor_FDR)]
  marker.table = as.data.table(read.table(marker.file, sep = "\t",header = T))
  marker.table = merge.data.table(p2geneDF,marker.table, by = "coordinates", all.x = TRUE)
  marker.table = marker.table[!is.na(peak.id) & !is.na(Correlation),]
  write.table(marker.table,file = marker.file, sep = "\t")
}
