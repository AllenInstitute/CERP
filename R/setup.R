#' Setup ArchR for genome and threads
#'
#' Sets global variables about specific genome build and number of threads to utilze on the machine.
#'
#' @param archr.genome Valid genome name for ArchR: mm10, hg38
#' @param archr.threads Number of threads to utilize, warning increases memory usage as well.
#' 
#' @return None
#'
#' @keywords internal
setup = function(archr.genome, archr.threads){
    ## Genome annotation info using ArchR defaults
    addArchRGenome(genome = archr.genome)
    ## Multi-core
    addArchRThreads(threads = archr.threads) 
}
