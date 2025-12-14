#' @import S7
NULL

#' Create a read-only S7 property
#'
#' @param name The name of the property (string). Must match the slot name.
#' @param type The S7 class definition (e.g., class_numeric, class_Date).
#' @param cast Optional function to transform input (e.g., as.Date).
#' 
#' @name md_report
#' @return A list of markdown content used for producing the pdf report. 
md_report <- new_generic("md_report", "x", function(x, ...) {
    S7_dispatch()
})

#' Evaluate Analysis Quality
#' 
#' A S7 generics to print out the analysis summary.
#' @name evaluate
#' @export
evaluate <- new_generic("evaluate", "x", function(x, ...) {
    S7_dispatch()
})

#' Report Analysis Result
#' 
#' A S7 generics to print out the analysis summary.
#' @name report
#' @export
report <- new_generic("report", "x", function(x, ...) {
    S7_dispatch()
})
