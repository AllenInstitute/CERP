#' Calls peaks currently utilizing ArchR
#'
#' Sets global variables about specific genome build and number of threads to utilze on the machine.
#'
#' @param archr.proj An already setup ArchR project or a name (string) to build a new ArchR project. If passing an already defined archr.proj then arrow.file.dir and cell.annotation.file will not be used.
#' @param archr.genome Valid genome name for ArchR: mm10, hg38
#' @param groupBy Metadata field to group cells by
#' @param dataset Information about the dataset that these marker peaks belong to. E.g. Brain region, etc.
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
#' @param publish Set to 'TRUE' if you are ready to make the marker peak table viewable on Shiny.
#' @param ucsc.user A ucsc genome browser user/account which to contains the session to build links from.
#' @param ucsc.session A registred session for the UCSC genome browser user/account specified in `ucsc.user`.
#' 
#' @return ArchRProject, markerPeaks
#'
#' @export
peakCaller = function(archr.proj, archr.genome, groupBy, dataset, archr.threads=4,               
                      arrow.file.dir=NULL, cell.annotation.file=NULL,
                      archr.clustering=FALSE, varFeatures=15000, resolution=c(0.2,1,2), 
                      tileSize=25, normMethod="ReadsInTSS", maxCells=NULL,              
                      archr.visualize=FALSE, output.folder=NULL, publish=FALSE, ucsc.user=NULL, ucsc.session=NULL,
                      counts_matrix = NULL, calculate_p2gene = FALSE, samp.dat = NULL,max_distance = 250000
){
  
  ## Error handling
  archr.proj = tryCatch({
    
    ## Run basic checks before getting to far.
    #.run_checks()
    
    ## unique_peak_calling_id
    start_time = Sys.time()
    unique.id = gsub(start_time, pattern = " |PDT|-|:", replacement = "")
    print(paste0("unique.id for this run:", unique.id))
    
    ## Remove old directory and create new directory under ArchR project for marker peak results
    print("creating marker_peaks directory")
    #unlink(file.path(getOutputDirectory(archr.proj), "MarkerPeaks", groupBy))
    dir.create(file.path(getOutputDirectory(archr.proj), "MarkerPeaks", groupBy), showWarnings = TRUE, recursive=TRUE)
    #dir.create(file.path(getOutputDirectory(archr.proj), "MarkerPeaks"), showWarnings = TRUE)
    #dir.create(file.path(getOutputDirectory(archr.proj), "MarkerPeaks", groupBy), showWarnings = TRUE)
    
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
    addArchRThreads(1)
    archr.proj = groupCoverage(archr.proj = archr.proj,
                               groupBy = groupBy)
    
    ##
    print("Calling peaks")
    addArchRThreads(archr.threads)
    archr.proj = peakCalling(archr.proj = archr.proj,
                             groupBy = groupBy)
    
    ##
    print("Identifying marker peaks")
    addArchRThreads(archr.threads)
    archr.proj = markerPeaks(archr.proj = archr.proj,
                             groupBy = groupBy, 
                             unique.id = unique.id)
    

    ## Record marker table location 
    archr.proj@projectMetadata[[paste0("loc_marker_table_", groupBy)]] = file.path(getOutputDirectory(archr.proj), "MarkerPeaks", groupBy, paste0(groupBy, "_markerPeaks.tsv"))
    
    ##
    if(!is.null(archr.proj@projectMetadata[[paste0(groupBy, "_markerPeaks")]])){
      print("Loading from project")
      marker.table = archr.proj@projectMetadata[[paste0(groupBy, "_markerPeaks")]]
    }else{
      print("Loading from csv")
      marker.table = read.csv(archr.proj@projectMetadata[[paste0("loc_marker_table_", groupBy)]])
    }
    
    ##
    annotatePeaks(marker.table = marker.table, 
                  archr.proj = archr.proj, 
                  groupBy = groupBy, 
                  dataset = dataset,
                  publish = publish,
                  filename = file.path(file.path(getOutputDirectory(archr.proj), "MarkerPeaks", groupBy), paste0(groupBy, "_annotated_markerPeaks.tsv")),
                  ucsc.user = ucsc.user,
                  ucsc.session = ucsc.session)


    if(calculate_p2gene == TRUE){
    if("GeneExpressionMatrix" %in% getAvailableMatrices(archr.proj) & 'LSI_Combined' %in% names(archr.proj@reducedDims)){
      print("Calculating peak to gene links")
      peak2Gene(archr.proj = archr.proj , max_distance = max_distance, reduced_dims = "LSI_Combined", useMatrix = "GeneExpressionMatrix",
                     marker.file = file.path(file.path(getOutputDirectory(archr.proj), "MarkerPeaks", groupBy), paste0(groupBy, "_annotated_markerPeaks.tsv")))
    }else{print("necessary ArchR elements not found to run peak to gene correlations")}

    }
    
    ##
    addArchRThreads(1)
    print("Producing bigwig and fragment files")
    generateBigWigs(archr.proj = archr.proj,
                    groupBy = groupBy)
    
    ## Return the archr.proj outside of try-catch
    archr.proj
  }, error=function(error) {
    message(error)
    return(archr.proj)
  }, finally={ 
    print("Done!") 
    ## Save project
    saveArchRProject(archr.proj)
  }
  )
  return(archr.proj) 
}
