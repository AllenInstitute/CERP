#' Check and sanatize user input
#'
#' @keywords int
.run_checks(){
    .checkMACS2() ## Check that we can load peak caller
}

#' Check that MACS2 can be found early on to save users time if MACS2 does not exist
#'
#' @keywords internal
.checkMACS2(){
    findMacs2()
}