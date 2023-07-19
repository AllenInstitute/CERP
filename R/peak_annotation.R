#' Post-pipeline annotation of marker peaks
#'
#' Defines additional marker peak annotations and standardizes reporting tables.
#'
#' @param marker.table file.path to a tab-seperated file containing results of ArchR marker peak analysis
#' @param archr.proj An already setup ArchR project or a string to build a new ArchR project.
#' @param dataset Information about the dataset that these marker peaks belong to. E.g. Brain region, etc.
#' @param taxonomy A defined taxonomy ID assocaiated with cell type annotations used for peak calling.
#' @param groupBy The annotation level used to group data for peak calling. E.g. Class, Subclass, etc.
#' @param publish Set to 'TRUE' if you are ready to make the marker peak table viewable on Shiny.
#' @param filename filename to save annotated marker table
#' @param ucsc.user A ucsc genome browser user/account which to contains the session to build links from.
#' @param ucsc.session A registred session for the UCSC genome browser user/account specified in `ucsc.user`.
#' 
#' @return annotated marker table
#'
#' @export
annotatePeaks = function(marker.table, archr.proj, dataset, taxonomy, groupBy, publish, filename, shiny.dir=NULL, ucsc.user=NULL, ucsc.session=NULL){

    ##
    marker.table$genome_coordinates = paste0(marker.table$seqnames, ":", marker.table$start, "-", marker.table$end)

    ## Get peak GRanges object from archr.project and order accorind to marker peaks
    peak.gr = getPeakSet(archr.proj) %>% as.data.frame(row.names=NULL)
    peak.gr$id = paste0(peak.gr$seqnames, ":", peak.gr$start, "-", peak.gr$end)
    peak.gr = peak.gr[match(marker.table$genome_coordinates, peak.gr$id),]

    ## Sanity check on ordering
    stopifnot(all(peak.gr$id == marker.table$genome_coordinates)) ##

    ## -- Annotate -------------

    ## Start reporting table with peak ID
    marker.reporting = data.frame(peak.id = paste0(getArchRGenome(), "_", dataset, "_", taxonomy, "_", groupBy, "_", archr.proj$unique.id, ":peak-", marker.table$idx))

    ## Cell population 
    marker.reporting$cell.population = marker.table$group_name

    ## 
    marker.reporting$coordinates = marker.table$genome_coordinates
    marker.reporting$chr = marker.table$seqnames
    marker.reporting$start = marker.table$start
    marker.reporting$end = marker.table$end
    marker.reporting$width = marker.table$end - marker.table$start
    marker.reporting$rank = marker.table %>% group_by(marker.table$group_name) %>% mutate(id = row_number()) %>% pull(id)
    marker.reporting$fdr = marker.table$FDR
    marker.reporting$fold.change = marker.table$Log2FC
    marker.reporting$peak.magnitude = peak.gr$score

    ## Gene info
    marker.reporting$nearest.gene = peak.gr$nearestGene
    marker.reporting$dist.to.gene = peak.gr$distToGeneStart
    marker.reporting$relative.position.gene = ""

    ## TSS info
    marker.reporting$nearest.TSS = peak.gr$nearestTSS
    marker.reporting$dist.to.TSS = peak.gr$distToTSS

    ## GC content, peak content
    marker.reporting$GC.content = peak.gr$GC
    marker.reporting$peak.type = peak.gr$peakType

    ## Add browser track link outs
    if(is.null(ucsc.user)){
        marker.reporting$genome_track = ""
    }else{
        marker.reporting$genome_track = paste0("https://genome.ucsc.edu/s/", ucsc.user, "/", ucsc.session, "?position=", marker.reporting$chr, ":", marker.reporting$start, "-", marker.reporting$end)
    }

    ##
    if(publish == TRUE){
        write.table(marker.reporting, row.names=FALSE, col.names=TRUE, file=paste0(dataset, "_", taxonomy, "_", groupBy, "_annotated_markerPeaks.tsv"))
    }

    ## Save the 
    if(!is.null(filename)){
        write.table(marker.reporting, row.names=F, col.names=T, sep="\t", file=filename)
    }

    ##
    if(is.null(archr.proj@projectMetadata$markerAnnotatedPeaks)){ archr.proj@projectMetadata$markerAnnotatedPeaks = list() }
    archr.proj@projectMetadata$markerAnnotatedPeaks[[paste0(dataset, "_", taxonomy, "_", groupBy, "_annotated_markerPeaks.tsv")]] = marker.reporting

    return(archr.proj)
}