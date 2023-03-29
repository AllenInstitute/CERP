#' Gets peak to gene links and merge with existing marker peaks
#'
#' @notes archr.proj must have ar_id and ar_directory for this to work properly. Future updates will allow for counts matrices to generate gene expression matrices.
#'
#' @param archr.proj 
#' @param max_distance distance from the gene to search for linked peaks - default 250kb
#' @param reduced_dims reduced dimension parameter for p2gene links ArchR function, default "LSI_Combined"
#' @param useMatrix' matrix to use for constructing p2gene links, can be GeneExpressionMatrix (default) or GeneIntegrationMatrix
#' @param marker.file' annotated marker peaks file path
#'
#' @keywords internal
peak_to_gene = function(archr.proj = NULL, max_distance = 250000, reduced_dims = "LSI_Combined", useMatrix = "GeneExpressionMatrix",
                        marker.file = file.path(getOutputDirectory(archr.proj), "MarkerPeaks", groupBy, paste0(groupBy, "_annotated_markerPeaks.tsv"))){

    ## Ensure annotated peak file exists
    if(!file.exists(marker.file)){print("No peak annotation file exists, please check that CERP has run properly"); return(archr.proj)}
    
    ## Check for GEX and LSI if not then setup.
    if(!("GeneExpressionMatrix" %in% getAvailableMatrices(archr.proj) & 'LSI_Combined' %in% names(archr.proj@reducedDims))){
        if("ar_id" %in% colnames(archr.proj@cellColData)& "ar_directory" %in% colnames(archr.proj@cellColData)){
            archr.proj = add_gex_and_LSI(archr.proj = archr.proj)  
        }else{
            print("The ArchR project is not properly configured for peak2gene. Please refere to documentation and/or examples at: https://github.com/AllenInstitute/CERP")
        }
    }

    ## Compute peak to gene links
    print("Calculating peak to gene links")
    archr.proj = addPeak2GeneLinks(ArchRProj = archr.proj, reducedDims = reduced_dims, useMatrix = useMatrix, maxDist = max_distance)

    ## Build peak to gene table
    archr.proj = make_p2g_table(archr.proj = archr.proj,
                                reducedDims = reduced_dims,
                                useMatrix = useMatrix,
                                maxDist = max_distance,
                                marker.file = marker.file) 

    return(archr.proj)
}

#' Construct peak to gene link table
#'
#' @param archr.proj
#' @param marker.file
#'
#' @keywords internal
make_p2g_table = function(archr.proj){
    ## Gather peak to gene information
    p2geneDF = metadata(archr.proj@peakSet)$Peak2GeneLinks
    p2geneDF$correlated_gene <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
    p2geneDF$coordinates <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), ":", start(.), "-", end(.))})[p2geneDF$idxATAC]
    p2geneDF$cor_FDR<-p2geneDF$FDR
    p2geneDF = as.data.table(p2geneDF)
    p2geneDF = p2geneDF[,.(coordinates,correlated_gene,Correlation,cor_FDR)]
    p2geneDF = na.omit(p2geneDF)
    ## Load in annotated marker peak table and update
    write.table(p2geneDF, file = file.path(getOutputDirectory(archr.proj), "MarkerPeaks", groupBy), paste0(groupBy, "_peak_to_gene.tsv"), sep = "\t")
}

#' Add gene expression and joint LSI to archr project
#'
#' @param archr.proj
#'
#' @keywords internal
add_gex_and_LSI = function(archr.proj){
    print("generating Gene Expression Matrix using AR directories")
    
    ## Add gene expression data
    h5paths = paste0(unique(archr.proj$ar_directory), unique(archr.proj$ar_id), "/", 'outs/filtered_feature_bc_matrix.h5')
    seRNA = import10xFeatureMatrix(input = h5paths, names = unique(archr.proj$Sample),strictMatch = F)
    seRNA = seRNA[,colnames(seRNA) %in% archr.proj$cellNames]
    archr.proj <- addGeneExpressionMatrix(input = archr.proj, seRNA = seRNA, force = TRUE, strictMatch = T)  ##verror unused argument strictMatch = F strictMatch

    # Perform LSI on both modalities and then combine -------------------------
    print("Calculating LSI_ATAC")
    archr.proj <- addIterativeLSI(
        ArchRProj = archr.proj, 
        clusterParams = list(
            resolution = c(0.2, 1, 2), 
            sampleCells = 10000,
            n.start = 10
        ),
        saveIterations = TRUE,
        useMatrix = "TileMatrix", 
        depthCol = "nFrags",
        varFeatures = 15000,
        name = "LSI_ATAC"
    )

    ##
    print("Calculating LSI_RNA")
    archr.proj <- addIterativeLSI(
        ArchRProj = archr.proj, 
        clusterParams = list(
            resolution = c(0.2, 1, 2), 
            sampleCells = 10000,
            n.start = 10
        ),
        saveIterations = TRUE,
        useMatrix = "GeneExpressionMatrix", 
        depthCol = "Gex_nUMI",
        varFeatures = 2500,
        firstSelection = "variable",
        binarize = FALSE,
        name = "LSI_RNA"
    )
    print("Combining LSIs")
    archr.proj <- addCombinedDims(archr.proj, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")
    return(archr.proj)
}