#' pheman
#'
#' Create Manhattan plots for GWAS
#' Note: There is an issue with dev.off() if using RStudio
#' Dependencies: ggplot2
#' Suggested: RColorBrewer, ggrepel
#' @param d data frame, if not plato or plink format, must contain PHE, SNP, CHR, POS, pvalue columns, optional Shape
#' @param format format of input
#' @param groups optional file for organizing phenotypes
#' @param line optional pvalue threshold to draw red line at
#' @param annotate_snp list of RSIDs to annotate
#' @param annotate_p pvalue threshold to annotate
#' @param title optional string for plot title
#' @param morecolors boolean, expand color palette, requires RColorBrewer package, default FALSE
#' @param file file name of saved image
#' @param hgt height of plot in inches
#' @param wi width of plot in inches
#' @param res resolution of plot in pixels per inch
#' @return png image
#' @export
#' @examples
#' pheman(d, format, line, groups, annotate_snp, annotate, title, morecolors, file, hgt, wi, res)

pheman <- function(d, format="plotman", groups=NULL, line, annotate_snp, annotate_p, title=NULL, morecolors=FALSE, file="gman", hgt=7, wi=12, res=300 ){
  if (!requireNamespace(c("ggplot2"), quietly = TRUE)==TRUE) {
    stop("Please install ggplot2 to create visualization.", call. = FALSE)
  } else {
    require("ggplot2", quietly = TRUE)
  }

  if(format=="plato"){
    stop("PLINK format coming soon...")
  } else if(format=="plato-codom"){
    stop("PLATO format coming soon...")
  }

  #Sort data
  d$CHR <- factor(d$CHR, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"))
  d_order <- d[order(d$CHR, d$POS), ]
  d_order$pos_index <- seq.int(nrow(d_order))

  d_order_sub <- d_order[colnames(d_order) %in% c("SNP", "CHR", "POS", "pvalue", "pos_index")]


}
