setMethod(
  f = "show",
  signature = "bsaR",
  definition = function(object){

    if( ncol(object@AD) > 0 ){
      nsamp <- ncol(object@AD)
    } else {
      nsamp <- 0
    }
    nchrom <- length(unique(object@meta$CHROM))
    nvar <- nrow(object@meta)

    cat("***** Object of Class bsa *****\n")
    cat( paste( nsamp, "samples\n") )
    cat( paste( nchrom, "CHROMs\n") )
    cat( paste( format(nvar, big.mark=","), "variants\n") )
    cat("*****        *****         *****\n")
  }
)

setMethod(f = "$",
          signature = "bsaR",
          definition = function(x, name){
            ## 'name' is a character(1)
            slot(x, name)
          })

