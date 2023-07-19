#' Calls peaks currently utilizing ArchR
#'
#' Sets global variables about specific genome build and number of threads to utilze on the machine.
#'
#' @param archr.proj An already setup ArchR project or a name (string) to build a new ArchR project. If passing an already defined archr.proj then arrow.file.dir and cell.annotation.file will not be used.
#' @param groupBy Metadata field to group cells by
#' @param dataset Information about the dataset that these marker peaks belong to. E.g. Brain region, etc.
#' @param taxonomy A defined taxonomy ID assocaiated with cell type annotations used for peak calling.
#' @param genomeSize Number of threads to utilize, warning increases memory usage as well.
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
#' @param calculate_gini_index ADD INFO
#' @param calculate_peak2gene ADD INFO
#' @param max_distance ADD INFO
#' 
#' @return ArchRProject, markerPeaks
#'
#' @export
peakCaller = function(archr.proj, groupBy, dataset, taxonomy, 
                      genomeSize = NULL,              
                      arrow.file.dir=NULL, 
                      cell.annotation.file=NULL,
                      ## Clustering parameters
                      archr.clustering=FALSE, 
                      varFeatures=15000, 
                      resolution=c(0.2,1,2), 
                      ## Normalization
                      tileSize=25, 
                      normMethod="ReadsInTSS", 
                      maxCells=NULL,              
                      archr.visualize=FALSE){

    ## Error handling
    archr.proj = tryCatch({

            ## Run basic checks before getting to far.
            #.run_checks()

            ## unique_peak_calling_id
            start_time = Sys.time()
            archr.proj$unique.id = gsub(start_time, pattern = " |PDT|-|:", replacement = "")
            print(paste0("unique.id for this run:", unique(archr.proj$unique.id)))

            ## Create new directory under ArchR project for marker peak results
            dir.create(file.path(getOutputDirectory(archr.proj), "MarkerPeaks", groupBy), showWarnings = TRUE, recursive=TRUE)

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
                                     groupBy = groupBy,
                                     genomeSize=genomeSize)
            saveArchRProject(archr.proj)

            ##
            print("Identifying marker peaks")
            archr.proj = markerPeaks(archr.proj = archr.proj,
                                     groupBy = groupBy)

            ## Record marker table location 
            archr.proj@projectMetadata[[paste0("loc_marker_table_", groupBy)]] = file.path(getOutputDirectory(archr.proj), "MarkerPeaks", groupBy, paste0(groupBy, "_markerPeaks.tsv"))
            saveArchRProject(archr.proj)

            ##
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

#' Annotate peaks using standard CERP format
#'
#' @param archr.proj An already setup ArchR project or a name (string) to build a new ArchR project. If passing an already defined archr.proj then arrow.file.dir and cell.annotation.file will not be used.
#' @param groupBy Metadata field to group cells by
#' @param dataset Information about the dataset that these marker peaks belong to. E.g. Brain region, etc.
#' @param publish Set to 'TRUE' if you are ready to make the marker peak table viewable on Shiny.
#' @param ucsc.user A ucsc genome browser user/account which to contains the session to build links from.
#' @param ucsc.session A registred session for the UCSC genome browser user/account specified in `ucsc.user`.
#' 
#' @return ArchRProject
#'
#' @export
peakAnnotation = function(archr.proj, groupBy, dataset, publish, ucsc.user, ucsc.session){
    ##
    if(!is.null(archr.proj@projectMetadata[[paste0(groupBy, "_markerPeaks")]])){
        print("Loading from project")
        marker.table = archr.proj@projectMetadata[[paste0(groupBy, "_markerPeaks")]]
    }else{
        print("Loading from tsv")
        marker.table = read.table(archr.proj@projectMetadata[[paste0("loc_marker_table_", groupBy)]], sep = "\t", header = T)
    }

    ##
    archr.proj = annotatePeaks(marker.table = marker.table, 
                                archr.proj = archr.proj, 
                                groupBy = groupBy, 
                                dataset = dataset,
                                publish = publish,
                                filename = file.path(file.path(getOutputDirectory(archr.proj), "MarkerPeaks", groupBy), paste0(groupBy, "_annotated_markerPeaks.tsv")),
                                ucsc.user = ucsc.user,
                                ucsc.session = ucsc.session)
    
    ##
    return(archr.proj)
}

#' Calls peaks currently utilizing ArchR
#'
#' Sets global variables about specific genome build and number of threads to utilze on the machine.
#'
#' @param archr.proj An already setup ArchR project.
#' @param groupBy Metadata field to group cells by
#' @param max_sample_size 
#' 
#' @return ArchRProject
#'
#' @export
peakSpecificity = function(archr.proj, groupBy, max_sample_size=500, numWorkers=16){
    ## Calculate gini index for each marker peak
    archr.proj = calculate_gini_index(archr.proj = archr.proj,
                                        groupBy = groupBy,
                                        filename = file.path(file.path(getOutputDirectory(archr.proj), "MarkerPeaks", groupBy), paste0(groupBy, "_annotated_markerPeaks.tsv")),
                                        numWorkers = numWorkers,
                                        max_sample_size = max_sample_size)
    return(archr.proj)
}

#' Rank peaks using the PeakRankR approach
#'
#' Adds new ranking columns to the peak annotation table
#'
#' @param archr.proj An ArchR project after running through CERP
#' @param groupBy Metadata field to group cells by
#' @param marker.peak.table The annotated marker peak table from `peakAnnotation()`
#' @param bw.table A table with bigwig file names and corresponding annotations that match cell.population in marker.peak.table.
#' 
#' @return ArchRProject
#'
#' @export
peakRankeR_annotation = function(archr.proj, groupBy, marker.peak.table, bw.table=NULL){

    ##
    if(is.null(bw.table)){
        bw.files = list.files(file.path(getOutputDirectory(archr.proj), "GroupBigWigs", groupBy), pattern = '*.bw', full.names=T)
        bw.names = unlist(lapply(strsplit(gsub(paste0(file.path(getOutputDirectory(archr.proj), "GroupBigWigs", groupBy), "/"), "", bw.files), "-"), "[[", 1))
        bw.table = data.frame(bw_path = bw.files, sample_id = bw.names)
    }

    ##
    ranked.marker.peak.table = Peak_RankeR(tsv_file_df = marker.peak.table,
                                            group_by_column_name = "cell.population",
                                            background_group     = unique(marker.peak.table$cell.population),
                                            bw_table             = bw.table, 
                                            rank_sum             = TRUE,
                                            weights              = c(1,1,1))

    ## Save marker peaks
    write.table(marker.list, 
                file=file.path(getOutputDirectory(archr.proj), "MarkerPeaks/", groupBy, "/", groupBy, "_peakRankR_markerPeaks.tsv"), 
                sep="\t", 
                row.names=FALSE)

    ##
    return(archr.proj)
}

#' Calls peaks currently utilizing ArchR
#'
#' Sets global variables about specific genome build and number of threads to utilze on the machine.
#'
#' @param archr.proj An already setup ArchR project.
#' @param max_distance ADD INFO
#' @param reduced_dims 
#' @param useMatrix
#'
#' @return ArchRProject
#'
#' @export
peakToGene = function(archr.proj, max_distance, reduced_dims="LSI_Combined", useMatrix="GeneExpressionMatrix"){
    ##
    archr.proj = peak_to_gene(archr.proj = archr.proj,
                                max_distance = max_distance, 
                                reduced_dims = reduced_dims, 
                                useMatrix = useMatrix)
    ## 
    return(archr.proj)
}