#### Class definitions. ###

#' @title bsaR class
#'
#' @slot meta data.frame
#' @slot AD matrix
#' @slot Depth matrix
#' @slot Freq matrix
#' @slot window list
#'
setClass(Class = "bsaR",
         slots = c(meta = "data.frame",
                   AD = "matrix",
                   Depth = "matrix",
                   Freq = "matrix",
                   Window = "list"))



