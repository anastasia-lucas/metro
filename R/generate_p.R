#' generate_p()
#'
#' Generate dummy data
#' @param group Boolean include Group column
#' @param shape Boolean include Shape column
#' @return png image(s)
#' @export
#' @examples
#' generate_p(color, shape)

generate_p <- function(color=FALSE, shape=FALSE){
  pvals <- data.frame(SNP=paste("rs", seq(1:5000), sep=""))
  pvals$CHR <- rep(c(1:22, "X", "Y"), length.out=5000, each=200)
  pvals$POS <- rep(seq(1, 10000, by = 200), length.out=5000)
  pvals$pvalue <- runif(n=5000)
  if(color==TRUE){
    pvals$Color <- rep(paste("Color", seq(1:6), sep="") , length.out=5000, each=200)
  }
  if(shape==TRUE){
    pvals$Shape <- rep(paste("S", seq(1:3), sep="") , length.out=5000, each=1)
  }
  return(pvals)
}
