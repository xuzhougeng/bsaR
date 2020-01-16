
#' Create bsa object from vcf
#'
#' @param vcf.file vcf file path
#'
#' @importFrom vcfR read.vcfR getCHROM getPOS extract.gt
#'
#' @export
CreateBsaFromVcf <- function(vcf.file){

  # create new BSA object
  bsa <- new("bsaR")

  vcf <- read.vcfR(vcf.file)

  meta <- data.frame(CHROM = getCHROM(vcf),
                     POS   = getPOS(vcf),
                     biallele = is.biallelic(vcf))
  slot(bsa, "meta") <- meta
  slot(bsa, "AD") <- extract.gt(vcf, element = "AD")

  return(bsa)

}

#' calculate the depth of each sample
#'
#' @param bsa bsa object
#'
#' @export
CalcDepth <- function(bsa){

  x <- bsa@AD

  AD_COUNT_LIST <-  strsplit(x, ",", fixed = TRUE, useBytes = TRUE)
  AD_COUNT <- as.integer(unlist(AD_COUNT_LIST))
  Ref_COUNT <- AD_COUNT[seq.int(1,length(AD_COUNT),2)]
  Alt_COUNT <- AD_COUNT[seq.int(2,length(AD_COUNT),2)]
  AD_Depth <- Ref_COUNT + Alt_COUNT
  AD_Depth <- matrix(AD_Depth, nrow = nrow(x))
  row.names(AD_Depth) <- row.names(x)
  colnames(AD_Depth) <- colnames(x)

  slot(bsa, "Depth") <- AD_Depth

  return(bsa)

}

#' calculate the allele frequency of ref
#'
#' @param bsa bsaR object
#'
#' @export
CalcAltFreq <- function(bsa){

  x <- bsa@AD

  AD_COUNT_LIST <-  strsplit(x, ",", fixed = TRUE, useBytes = TRUE)
  AD_COUNT <- as.integer(unlist(AD_COUNT_LIST))
  Ref_COUNT <- AD_COUNT[seq.int(1,length(AD_COUNT),2)]
  Alt_COUNT <- AD_COUNT[seq.int(2,length(AD_COUNT),2)]
  AD_freq <- Alt_COUNT / (Ref_COUNT + Alt_COUNT)
  AD_freq <- matrix(AD_freq, nrow = nrow(x))
  row.names(AD_freq) <- row.names(x)
  colnames(AD_freq) <- colnames(x)

  slot(bsa, "Freq") <- AD_freq
  return(bsa)

}

#' caculate the frequency by windows
#'
#' @param bsa bsaR object
#' @param window.size window size, default is 50000
#'
#' @export
CalcFreqByWindow <- function(bsa, window.size = 50000){

  list.name <- paste0("w", window.size)

  if ( list.name %in% names(slot(bsa, "Window"))){
    return(bsa)
  }

  Freq <- bsa$Freq
  meta.split <- split(bsa$meta, bsa$meta$CHROM)
  variants.num <- vapply(meta.split, nrow,1)
  variants.cumsum <- cumsum(variants.num)

  freq.split <- lapply(unique(bsa$meta$CHROM), function(x) {

    pos <- meta.split[[x]]$POS
    ranges <- range(pos)
    # too short
    if ( ranges[2] - ranges[1] < 3 * window.size ) return(NULL)

    # break windows
    windows <- cut(pos, breaks = seq(1, ranges[2] + window.size, window.size))
    windows.levels <- levels(windows)

    # get the start and end of each chromosome
    start <- variants.cumsum[x] - variants.num[x] + 1
    end   <- variants.cumsum[x]
    # subset the Freq
    freqs <- Freq[seq(start, end), ]

    freq.windows <- matrix(0,nrow = length(windows.levels),ncol = 2)
    row.names(freq.windows) <- windows.levels
    colnames(freq.windows) <- colnames(Freq)

    for (window in windows.levels) {
      idx <- which(windows == window)
      freq.windows[window,] <- colMeans(freqs[idx, , drop = FALSE])
    }
    freq.windows
  })
  names(freq.split) <- unique(bsa$meta$CHROM)

  slot(bsa, "Window")[[list.name]] <- freq.split
  return(bsa)

}
