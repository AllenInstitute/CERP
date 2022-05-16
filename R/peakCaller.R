#' Calls peaks currently utilizing ArchR
#'
#' Sets global variables about specific genome build and number of threads to utilze on the machine.
#'
#' @param archr.proj An already setup ArchR project or a name (string) to build a new ArchR project. If passing an already defined archr.proj then arrow.file.dir and cell.annotation.file will not be used.
#' @param archr.genome Valid genome name for ArchR: mm10, hg38
#' @param groupBy Metadata field to group cells by
#' @param archr.threads Number of threads to utilize, warning increases memory usage as well.
#' @param arrow.file.dir Directory where all the Arrow files are stored 
#' @param cell.annotation.file File that contains sample_id and additional cell annotations to append to the ArchRProject
#' @param archr.clustering Boolean (T/F) determining if clustering/LSI should be run.
#' @param varFeatures Number of variable features to use for cluster
#' @param resolution Vector of clustering resolutions for multi-round clusters
#' @param tileSize Width of tiles that will section the genome
#' @param normMethod Normalization method for fragment pileups 
#' @param maxCells Number of cells to use for each group.
#' @param archr.visualize Should marker peak visualizations be produced
#' @param output.folder If supplied we will save the marker peak table and ArchRProject to this location.
#' 
#' @return ArchRProject, markerPeaks
#'
#' @export
peakCaller = function(archr.proj, archr.genome, groupBy, archr.threads=4,               
                      arrow.file.dir=NULL, cell.annotation.file=NULL,
                      archr.clustering=FALSE, varFeatures=15000, resolution=c(0.2,1,2), 
                      tileSize=25, normMethod="ReadsInTSS", maxCells=NULL,              
                      archr.visualize=TRUE, output.folder=NULL                                
                      ){

    ## Run basic checks before getting to far.
    .run_checks()

    ##
    print("Running initial ArchR setup")
    setup(archr.genome = archr.genome,
          archr.threads = archr.threads)

    ##
    print("Building/checking ArchRProject")
    archr.proj = build(arrow.file.dir = arrow.file.dir, 
                        cell.annotation.file = cell.annotation.file, 
                        archr.proj = archr.proj)

    ##
    if(archr.clustering == TRUE){
        print("Performing clustering and LSI")
        archr.proj = clusterDimReduc(archr.proj = archr.proj, 
                                    varFeatures = varFeatures, 
                                    resolution = resolution)
    }

    ##
    print("Defining group coverages")
    archr.proj = groupCoverage(archr.proj = archr.proj,
                               groupBy = groupBy)

    ##
    print("Calling peaks")
    archr.proj = peakCalling(archr.proj = archr.proj,
                             groupBy = groupBy)

    ##
    print("Identifying marker peaks")
    markerPeaks(archr.proj = archr.proj,
                groupBy = groupBy)

    ##
    print("Producing bigwig and fragment files")
    generateBigWigs(archr.proj = archr.proj,
                    groupBy = groupBy)

    ##
    return(archr.proj)
}