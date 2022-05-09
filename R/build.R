#' Builds the primary archr project object
#'
#' Sets global variables about specific genome build and number of threads to utilze on the machine.
#'
#' @param arrow.file.dir Valid genome name for ArchR: mm10, hg38
#' @param cell.annotation.file Number of threads to utilize, warning increases memory usage as well.
#' @param archr.proj An already setup ArchR project or a string to build a new ArchR project.
#' @param ... Additional parameters to ArchRProject()
#'
#' @return None
#'
#' @export
build = function(arrow.file.dir, cell.annotation.file, archr.proj="peakCallingPipeline", ...){

    ## Build ArchRProject 
    if(is.string(archr.proj)){
        print("Building ArchR Project")
        
        ## Identify arrow files
        arrow.files = list.files(arrow.file.dir, pattern="*.arrow")

        ## Sanity check arrow files
        if(length(arrow.files) == 0){stop("No arrow files found, please check input.")}

        ## Assemble ArchR project from Arrow files
        archr.proj = ArchR::ArchRProject(
            ArrowFiles = arrow.files, 
            outputDirectory = archr.proj,
            copyArrows = TRUE,
            ...
        )
    }else{
        if(is.null(archr.proj)){ stop("archr.proj must not be NULL. Either pre-build ArchRProject or string for building a new ArchRProject")}
        
        ## Passed basic checks, using user supplied ArchRProject
        print("Using predefined ArchR Project")
    }

    ## Add metadata to ArchR project
    if(grep("csv", cell.annotation.file)){
        cell.annotations = read.csv(cell.annotation.file)
    }else if(grep("feather", cell.annotation.file)){
        cell.annotations = read.feather(cell.annotation.file)
    }else if(is.null(cell.annotation.file)){
        cell.annotations = NULL
    }else{
        stop("Unidentified file type for cell annotations, should be csv or feather format")
    }

    ## Add metadata 
    archr.proj = addCellColData(ArchRProj = archr.proj,
                                data = cell.annotations)

    ##
    return(archr.proj)
}