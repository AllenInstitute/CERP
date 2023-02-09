peakCaller = function(archr.proj, archr.genome, groupBy, dataset, archr.threads=4,               
                      arrow.file.dir=NULL, cell.annotation.file=NULL,
                      archr.clustering=FALSE, varFeatures=15000, resolution=c(0.2,1,2), 
                      tileSize=25, normMethod="ReadsInTSS", maxCells=NULL,              
                      archr.visualize=FALSE, output.folder=NULL, publish=FALSE, ucsc.user=NULL, ucsc.session=NULL,
                      calculate_p2gene = FALSE, samp.dat = NULL,max_distance = 250000
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
    
