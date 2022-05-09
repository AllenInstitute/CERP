#' Calls peaks currently utilizing ArchR
#'
#' Sets global variables about specific genome build and number of threads to utilze on the machine.
#'
#' @param arrow.file.dir Valid genome name for ArchR: mm10, hg38
#' @param cell.annotation.file Number of threads to utilize, warning increases memory usage as well.
#' @param archr.proj An already setup ArchR project or a string to build a new ArchR project.
#' @param groupBy Metadata field to group cells by
#' @param tileSize Width of tiles that will section the genome
#' @param normMethod Normalization method for fragment pileups 
#' @param maxCells Number of cells to use for each group.
#' @param varFeatures Number of variable features to use for cluster
#' @param resolution Vector of clustering resolutions for multi-round clusters
#' @param output.folder If specified, will save marker heatmap pdf to output folder location, along with the peak matrix
#'
#' @return ArchRProject
#'
#' @export
peakCaller = function(){

    ##
    setup()

    ##
    archr.proj = build()

    ##
    archr.proj = clusterDimReduc()

    ##
    archr.proj = groupCoverage()

    ##
    archr.proj = peakCalling()

    ##
    list = markerPeaks()

    ##
    generateBigWigs()

    return(archr.proj)
}