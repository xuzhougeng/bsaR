#' filter multi-variant loci
#'
#' @param bsa bsa object
#'
#' @export
FilterMultiVariant <- function(bsa){

  message("only keep bi-variant loci")
  mask <- bsa@meta$biallele

  slot(bsa, "meta") <- bsa@meta[mask, ]
  slot(bsa, "AD") <- bsa@AD[mask,]

  if (! nrow(bsa@Freq) == 0 ){
    slot(bsa, "Freq") <- bsa@Freq[mask,]
  }

  if (! nrow(bsa@Depth) == 0 ){
    slot(bsa, "Depth") <- bsa@Depth[mask,]
  }

  return(bsa)


}


#' filter NA Variant
#'
#' @param bsa bsa object
#'
#' @export
FilterNaVariant <- function(bsa){

  message("filter NA position")

  mask <- rowSums(is.na(bsa@AD)) == 0

  slot(bsa, "meta") <- bsa@meta[mask, ]
  slot(bsa, "AD") <- bsa@AD[mask,]

  if (! nrow(bsa@Freq) == 0 ){
    slot(bsa, "Freq") <- bsa@Freq[mask,]
  }

  if (! nrow(bsa@Depth) == 0 ){
    slot(bsa, "Depth") <- bsa@Depth[mask,]
  }

  return(bsa)


}


#' filter allele matrix by depth
#'
#' @param bsa bsa object
#' @param depth depth threshold, default is 20
#'
#' @export
FilterLowDepth <- function(bsa, depth = 20){

  mask <- apply(bsa@Depth >= depth, 1, all)

  slot(bsa, "meta") <- bsa@meta[mask, ]
  slot(bsa, "AD") <- bsa@AD[mask,]

  if (! nrow(bsa@Freq) == 0 ){
    slot(bsa, "Freq") <- bsa@Freq[mask,]
  }

  if (! nrow(bsa@Depth) == 0 ){
    slot(bsa, "Depth") <- bsa@Depth[mask,]
  }


  return(bsa)

}

#' subset genome
#'
#' @param chr chromosome
#' @param start start position
#' @param end end position
#'
#' @export
SubsetGenome <- function(bsa, chr, start = -Inf, end = Inf){

  if (end < start) stop(" end should be great than start ")

  mask <- bsa@meta$CHROM == chr &
    bsa@meta$POS >= start &
    bsa@meta$POS <= end

  slot(bsa, "meta") <- bsa@meta[mask, ]
  slot(bsa, "AD") <- bsa@AD[mask,]

  if (! nrow(bsa@Freq) == 0 ){
    slot(bsa, "Freq") <- bsa@Freq[mask,]
  }

  if (! nrow(bsa@Depth) == 0 ){
    slot(bsa, "Depth") <- bsa@Depth[mask,]
  }


  return(bsa)

}
