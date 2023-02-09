#'gets peak to gene links and merge with existing marker peaks'
# archr.proj must have ar_id and ar_directory for this to work properly. Future updates will allow for counts matrices to generate gene expression matrices.
#' @param archr.proj 
#' @param max_distance distance from the gene to search for linked peaks - default 250kb
#' @param reduced_dims reduced dimension parameter for p2gene links ArchR function, default "LSI_Combined"
#' @param useMatrix' matrix to use for constructing p2gene links, can be GeneExpressionMatrix (default) or GeneIntegrationMatrix
#' @param marker.file' annotated marker peaks file path
#' @param reduced_dims reduced dimension for p2gene links
#' @param useMatrix matrix to use for constructing p2gene links
#' @param maxDist maximum distance from TSS to look for correlated peaks


peak_to_Gene = function(archr.proj = NULL, max_distance = 250000, reduced_dims = "LSI_Combined", useMatrix = "GeneExpressionMatrix",
                     marker.file = file.path(file.path(getOutputDirectory(archr.proj), "MarkerPeaks", groupBy), paste0(groupBy, "_annotated_markerPeaks.tsv"))){
  #check for LSI if not there, run LSI 
  if("GeneExpressionMatrix" %in% getAvailableMatrices(archr.proj) & 'LSI_Combined' %in% names(archr.proj@reducedDims)){
    archr.proj = make_p2g_table(archr.proj = archr.proj,reducedDims = reduced_dims,useMatrix = useMatrix,maxDist = max_distance,marker.file = marker.file)
  }else if(("GeneExpressionMatrix" %in% getAvailableMatrices(archr.proj) & 'LSI_Combined' %in% names(archr.proj@reducedDims)) == FALSE){
    if("ar_id" %in% colnames(archr.proj@cellColData)& "ar_directory" %in% colnames(archr.proj@cellColData)){
      archr.proj = add_gex_and_LSI(archr.proj = archr.proj)
      archr.proj = make_p2g_table(archr.proj = archr.proj,reducedDims = reduced_dims,useMatrix = useMatrix,maxDist = max_distance,marker.file = marker.file)
  }
  }
  return(archr.proj)
}


make_p2g_table = function(archr.proj = NULL,reducedDims = NULL,useMatrix = NULL,maxDist = NULL,marker.file = NULL){
    print("calculating peak to gene links")
    archr.proj <- addPeak2GeneLinks(ArchRProj = archr.proj,reducedDims = reduced_dims,useMatrix = useMatrix,maxDist = max_distance)
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
    return(archr.proj)
}


add_gex_and_LSI = function(archr.proj = NULL ){
      print("generating Gene Expression Matrix using AR directories")
      h5paths = paste0(unique(archr.proj$ar_directory), unique(archr.proj$ar_id), "/", 'outs/filtered_feature_bc_matrix.h5')
      seRNA = import10xFeatureMatrix(input = h5paths, names = unique(archr.proj$Sample),strictMatch = F)
      seRNA = seRNA[,colnames(seRNA) %in% archr.proj$cellNames]
      archr.proj <- addGeneExpressionMatrix(input = archr.proj, seRNA = seRNA, force = TRUE,strictMatch = F)  #error unused argument strictMatch = F strictMatch
      
      # Perform LSI on both modalities and then combine -------------------------
      archr.proj <- addIterativeLSI(ArchRProj = archr.proj,
                                    clusterParams = list(resolution = c(0.2, 1, 2),
                                    sampleCells = 10000,n.start = 10),saveIterations = FALSE,
                                    useMatrix = "TileMatrix",depthCol = "nFrags",varFeatures = 15000,name = "LSI_ATAC")
      
      archr.proj <- addIterativeLSI(ArchRProj = archr.proj,
                                    clusterParams = list(resolution = c(0.2, 1, 2),
                                    sampleCells = 10000,n.start = 10),
                                    saveIterations = FALSE,
                                    useMatrix = "GeneExpressionMatrix",
                                    depthCol = "Gex_nUMI",
                                    varFeatures = 2500,
                                    firstSelection = "variable",
                                    binarize = FALSE,
                                    name = "LSI_RNA"
      )
      archr.proj <- addCombinedDims(archr.proj, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")
      return(archr.proj)
}

