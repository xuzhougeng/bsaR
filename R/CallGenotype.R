#' call geneotyep by allele frequency
#' 
#' The function will call genotype by allele frequency
#' 0 means homozygous reference, 1 means homozygous alternative, 2 means heterozygous
#' 
#' @param bsa bsaR object
#' @param low.AF low allele frequency threshold, default is 0.2
#' @param high.AF high allele frequency threshold, default is 0.8
#' 
#' @export 
CallGenotype <- function(bsa, low.AF=0.2, high.AF=0.8){
    x <- bsa@Freq
    # x should not be null
    if (length(x) == 0){
        stop("The allele frequency matrix is empty, run CalcAltFreq first")
    }
    # call genotype by allele frequency
    # if AF < low.AF, call 0, if AF > high.AF, call 1, otherwise call 2
    homo_ref_pos <- which(x < low.AF)
    homo_alt_pos <- which(x > high.AF)
    heter_pos <- which( x >= low.AF & x <= high.AF)
    x[homo_ref_pos] <- 0
    x[homo_alt_pos] <- 1
    x[heter_pos] <- 2

    slot(bsa, "GT") <- x

    return(bsa)
}