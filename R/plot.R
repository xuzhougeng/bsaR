#' plot window frequency
#'
#' @param bsa bsa object
#' @param geno1 genotype 1
#' @param geno2 genotype 2
#' @param chr chromosome
#' @param window.size widows.size
#'
#' @export
plotWindowFreq <- function(bsa,
                           chr = NULL,
                           geno1 = NULL,
                           geno2 = NULL,
                           window.size = NULL,
                           cols = c("red","blue","black")){
  list.name <- paste0("w", window.size)
  windows <- slot(bsa, "Window")[[list.name]]
  windows.names <- names(windows)

  if (is.null(chr)){
    stop("please select the contig/chromosome")
  }
  if ( ! chr %in% windows.names ){
    stop("The contig/chromosome name is wrong")
  }
  if ( is.null(geno1) ){
    stop("The geno1 is not set")
  }
  if ( is.null(window.size) ){
    stop("The window.size is not set")
  }

  chr_mt <- windows[[chr]]

  pos <- chr_mt[,"pos"]

  plot(x = pos, y=chr_mt[, geno1],
       xlab = "Position",
       ylab = "Delta",
       col = cols[1],
       pch = 19,
       cex = 0.5,
       ylim = c(-1,1),
       axes = FALSE,
       frame.plot = TRUE)
  if (!is.null(geno2)){
    points(x =pos, y=chr_mt[, geno2],
           col= cols[2], pch = 19, cex = 0.5)
  }
  delta <- chr_mt[, geno2] - chr_mt[, geno1]
  lines(x = pos, y = delta, col = cols[3])

  break_points <- seq(1, max(pos) + 1000000, 2000000)
  break_labels <- paste0(floor(break_points/1000000), "M")

  axis(1, at = break_points, labels = break_labels)
  axis(2, at = seq(-1, 1, 0.2))

}
