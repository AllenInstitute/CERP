#'gets peak to gene links and merge with existing marker peaks'

#' @param archr.proj 
#' @param max_distance distance from the gene to search for linked peaks - default 250kb
#' @param reduced_dims reduced dimension parameter for p2gene links ArchR function, default "LSI_Combined"
#' @param useMatrix' matrix to use for constructing p2gene links, can be GeneExpressionMatrix (default) or GeneIntegrationMatrix
#' @param marker.file' annotated marker peaks file path

