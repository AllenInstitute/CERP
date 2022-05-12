#' Check and sanatize user input
#'
#' @keywords internal
.run_checks = function(){
    .checkMACS2() ## Check that we can load peak caller
}

#' Check that MACS2 can be found early on to save users time if MACS2 does not exist
#'
#' @keywords internal
.checkMACS2 = function(){
    findMacs2()
}