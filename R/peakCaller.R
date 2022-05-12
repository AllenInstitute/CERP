#' Calls peaks currently utilizing ArchR
#'
#' Sets global variables about specific genome build and number of threads to utilze on the machine.
#'
#' @param archr.genome Valid genome name for ArchR: mm10, hg38
#' @param archr.threads Number of threads to utilize, warning increases memory usage as well.
#' @param arrow.file.dir Directory where all the Arrow files are stored 
#' @param cell.annotation.file File that contains sample_id and additional cell annotations to append to the ArchRProject
#' @param archr.proj An already setup ArchR project or a name (string) to build a new ArchR project. If passing an already defined archr.proj then arrow.file.dir and cell.annotation.file will not be used.
#' @param archr.clustering Boolean (T/F) determining if clustering/LSI should be run.
#' @param varFeatures Number of variable features to use for cluster
#' @param resolution Vector of clustering resolutions for multi-round clusters
#' @param groupBy Metadata field to group cells by
#' @param tileSize Width of tiles that will section the genome
#' @param normMethod Normalization method for fragment pileups 
#' @param maxCells Number of cells to use for each group.
#' @param output.folder If supplied we will save the marker peak table and ArchRProject to this location.
#'
#' @return ArchRProject, markerPeaks
#'
#' @export
peakCaller = function(archr.genome, archr.threads=4,                                    ## Setup
                      arrow.file.dir=NULL, cell.annotation.file=NULL, archr.proj,       ## Build
                      archr.clustering=FALSE, varFeatures=15000, resolution=c(0.2,1,2), ## clusterDimReduc
                      groupBy,                                                          ## groupCoverage, peakCalling, markerPeaks
                      tileSize=25, normMethod="ReadsInTSS", maxCells=NULL,              ## generateBigWigs
                      output.folder=NULL                                
                      ){

    ##
    setup(genome = archr.genome,
            archr.threads = archr.threads)

    ##
    archr.proj = build(arrow.file.dir, 
                        cell.annotation.file, 
                        archr.proj="peakCallingPipeline")

    ##
    if(archr.clustering == TRUE){
        archr.proj = clusterDimReduc(archr.proj, 
                                    varFeatures=varFeatures, 
                                    resolution=resolution)
    }
    ##
    archr.proj = groupCoverage(archr.proj,
                               groupBy)

    ##
    archr.proj = peakCalling(archr.proj,
                             groupBy)

    ##
    marker.peaks = markerPeaks(archr.proj,
                               groupBy)

    ##
    generateBigWigs(archr.proj,
                    groupBy)

    ##
    return.obj = list(archr.proj, marker.peaks); names(return.obj) = c("ArchRProject", "markerPeaks")
    return(return.obj)
}