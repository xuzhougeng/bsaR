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
#' @param min.depth minimum depth threshold, default is 20
#' @param max.depth maximum depth threshold, default is 3 * min.depth
#'
#' @export
FilterByDepth <- function(bsa, min.depth = 20, max.depth = 3 * min.depth){

  mask <- apply(bsa@Depth >= min.depth & bsa@Depth <= max.depth, 1, all)

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

#' filter allele matrix by low AF
#'
#' @param bsa bsa object
#' @param min.AF minimum AF threshold, default is 0.01
#' @param max.AF maximum AF threshold, default is 0.01
#'
#' @export
FilterByAF <- function(bsa, min.AF = 0.01, max.AF = 1 ){

  mask_mt <- bsa@Freq >= min.AF & bsa@Freq <= max.AF
  mask <- rowSums(mask_mt) == 2

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

#' filter loci  by geno type
#' 
#' @param bsa bsa object
#' @param p.wt sample name of parent wild type
#' @param p.mut sample name of parent mutant
#' 
#' @export 
FilterByGeno <- function(bsa, p.wt, p.mut){
  x <- bsa@GT
  # slot GT should be empty
  if (length(x) == 0){
    stop("GT slot should be empty, run CallGenotype first")
  }
  # the sample name should be in the GT matrix
  if (! p.wt %in% colnames(x) | ! p.mut %in% colnames(x)){
    stop("The sample name is not in the AD matrix, check your input")
  }
  #  the genotype of parent should not same
  p.wt.gt <- x[, p.wt]
  p.mut.gt <- x[, p.mut]
  mask1 <- p.wt.gt != p.mut.gt
  
  # the genotype of parent should not be 2
  mask2 <- p.wt.gt != 2 & p.mut.gt != 2

  mask <- which( mask1 & mask2)

  slot(bsa, "meta") <- bsa@meta[mask, ]
  slot(bsa, "Depth") <- bsa@Depth[mask,]
  slot(bsa, "AD") <- bsa@AD[mask,]
  slot(bsa, "GT") <- bsa@GT[mask,]
  slot(bsa, "Freq") <- bsa@Freq[mask,]

  return (bsa)

}


#' filter allele matrix by Position
#'
#' @param bsa bsa object
#' @param pos position vector, e.g. chr1_1000
#'
#' @export
FilterByPos <- function(bsa, pos ){

  all.pos <- paste0(bsa@meta$CHROM, "_", bsa@meta$POS)
  mask <- all.pos %in% conf_marker

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
