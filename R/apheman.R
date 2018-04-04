#' agman
#'
#' Create Manhattan plots for PheWAS
#' Note: There is an issue with dev.off() if using RStudio
#' Dependencies: ggplot2, gganimate
#' Suggested: RColorBrewer
#' @param d data frame, if not plato or plink format, must contain PHE, SNP, CHR, POS, pvalue, Frame columns, optional Shape
#' @param format format of input
#' @param phegroup optional grouping file for phenotypes, must contain PHE and Group columns
#' @param line optional pvalue threshold to draw red line at
#' @param log10 plot -log10() of pvalue column, boolean
#' @param yaxis label for y-axis, automatically set if log10=TRUE
#' @param opacity opacity of points, from 0 to 1, useful for dense plots
#' @param title optional string for plot title
#' @param chrcolor1 first alternating color for chromosome
#' @param chrcolor2 second alternating color for chromosome
#' @param highlight_snp list of snps to highlight
#' @param highlight_p pvalue threshold to highlight
#' @param highlighter color to highlight
#' @param groupcolors named list of colors for data in 'Color' column
#' @param file file name of saved image
#' @param ext file type to save, "gif" or "html"
#' @param hgt height of plot in inches
#' @param wi width of plot in inches
#' @param res resolution of plot in pixels per inch
#' @return png image
#' @export
#' @examples
#' apheman(d, format, phegroup, line, log10, yaxis, opacity, title, chrcolor1, chrcolor2, groupcolors, file, ext, hgt, wi, res)

apheman <- function(d, format="plotman", phegroup, line, log10=TRUE, yaxis, opacity=1, highlight_snp, highlight_p, highlighter="red", title=NULL, chrcolor1="#AAAAAA", chrcolor2="#4D4D4D", groupcolors, file="apheman", ext="html", hgt=7, wi=12, res=300 ){
  if (!requireNamespace(c("ggplot2"), quietly = TRUE)==TRUE|!requireNamespace(c("gganimate"), quietly = TRUE)==TRUE) {
    stop("Please install ggplot2 and ggiraph to create interactive visualization.", call. = FALSE)
  } else {
    require("ggplot2", quietly=TRUE)
    require("gganimate", quietly=TRUE)
  }

  #Check that frame is present
  if(!("Frame" %in% names(d))){stop("Please add 'Frame' column for animation attribute")}

  #Format input
  if(format=="plink"){
    stop("PLINK format coming soon...")
  } else if(format=="plato"){
    stop("PLATO format coming soon...")
  } else if(format=="plato-codom"){
    stop("PLATO-codom format coming soon...")
  }

  #Sort data
  d$CHR <- factor(d$CHR, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"))
  if(!missing(phegroup)){
    print("Only phenotypes with grouping information will be plotted")
    d_phe <- merge(phegroup, d, by="PHE")
    names(d_phe)[names(d_phe)=="Group"] <- "Color"
  } else {
    names(d_phe)[names(d_phe)=="PHE"] <- "Color"
  }
  d_order <- d_phe[order(d_phe$CHR, d_phe$POS), ]
  d_order$pos_index <- seq.int(nrow(d_order))
  d_order_sub <- d_order[colnames(d_order) %in% c("SNP", "CHR", "POS", "pvalue", "pos_index")]

  #Set up dataframe with color and position info
  maxRows <- by(d_order_sub, d_order_sub$CHR, function(x) x[which.max(x$pos_index),])
  minRows <- by(d_order_sub, d_order_sub$CHR, function(x) x[which.min(x$pos_index),])
  milimits <- do.call(rbind, minRows)
  malimits <- do.call(rbind, maxRows)
  lims <- merge(milimits, malimits, by="CHR")
  names(lims) <- c("Color", "snpx", "px", "posx", "posmin", "snpy", "py", "posy", "posmax")
  lims$av <- (lims$posmin + lims$posmax)/2
  lims <- lims[order(lims$Color),]
  lims$shademap <- rep(c("shade_ffffff", "shade_ebebeb"), each=1)

  #Frame must be a factor
  d_order$Frame <- factor(d_order$Frame)

  #Set up colors
  nchrcolors <- nlevels(factor(lims$Color))
  if(!missing(groupcolors)){
     dcols <- c(rep(x=c(chrcolor1, chrcolor2), length.out=nchrcolors, each=1), "#FFFFFF", "#EBEBEB")
     names(dcols) <-c(levels(factor(lims$Color)), "shade_ffffff", "shade_ebebeb")
    newcols <- c(dcols, groupcolors)
  } else {
    if (!requireNamespace(c("RColorBrewer"), quietly = TRUE)==TRUE) {
      stop("Please install RColorBrewer to add color attribute.", call. = FALSE)
    } else {
       require("RColorBrewer", quietly=TRUE)
    }
    ngroupcolors <- nlevels(factor(d_order$Color))
    getPalette = colorRampPalette(brewer.pal(11, "Spectral"))
    newcols <- c(rep(x=c(chrcolor1, chrcolor2), length.out=nchrcolors, each=1), getPalette(ngroupcolors), "#FFFFFF", "#EBEBEB")
    names(newcols) <-c(levels(factor(lims$Color)), levels(factor(d_order$Color)), "shade_ffffff", "shade_ebebeb")
  }

  #Info for y-axis
  if(log10==TRUE){
    d_order$pval <- -log10(d_order$pvalue)
    yaxislab <- expression(paste("-log"[10], "(p-value)", sep=""))
    if(!missing(line)) {redline <- -log10(line)}
  } else {
    d_order$pval <- d_order$pvalue
    yaxislab <- yaxis
    if(!missing(line)) {redline <- line}
  }

  #Start plotting
  p <- ggplot() + geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = 0, ymax = Inf, fill=factor(shademap)), alpha = 0.5)
  #Add shape info if available
  if("Shape" %in% names(d)){
    p <- p + geom_point(data=d_order, aes(x=pos_index, y=pval, color=factor(Color), shape=factor(Shape), frame=Frame), alpha=opacity)
  } else {
    p <- p + geom_point(data=d_order, aes(x=pos_index, y=pval, color=factor(Color), frame=Frame), alpha=opacity)
  }
  #if(!missing(annotate_p)){
  #  if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE) {
  #    print("Consider installing 'ggrepel' for improved text annotation")
  #    p <- p + geom_text(data=d_order[d_order$pvalue < annotate_p,], aes(pos_index,pval,label=SNP))
  #  } else {
  #    require("ggrepel", quietly = TRUE)
  #    p <- p + geom_text_repel(data=d_order[d_order$pvalue < annotate_p,], aes(pos_index,pval,label=SNP))
  #  }
  #}
  #if(!missing(annotate_snp)){
  #  if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE){
  #    print("Consider installing 'ggrepel' for improved text annotation")
  #    p <- p + geom_text(data=d_order[d_order$SNP %in% annotate_snp,], aes(pos_index,pval,label=SNP))
  #  } else {
  #    require("ggrepel", quietly = TRUE)
  #    p <- p + geom_text_repel(data=d_order[d_order$SNP %in% annotate_snp,], aes(pos_index,pval,label=SNP))
  #  }
  #}
  p <- p + scale_x_continuous(breaks=lims$av, labels=lims$Color, expand=c(0,0))
  p <- p + geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = -Inf, ymax = 0, fill=as.factor(Color)), alpha = 1)
  #Add legend
  p <- p + scale_colour_manual(name = "Color", values = newcols) + scale_fill_manual(name = "Color", values = newcols, guides(alpha=FALSE))
  p <- p + theme(axis.text.x=element_text(angle=90), panel.grid.minor.x = element_blank(), panel.grid.major.x=element_blank(), axis.title.x=element_blank(), legend.position="bottom", legend.title=element_blank())
  #Highlight if given
  if(!missing(highlight_snp)){
    if("Shape" %in% names(d)){
      p <- p + geom_point(data=d_order[d_order$SNP %in% highlight_snp, ], aes(x=pos_index, y=pval, shape=Shape), colour=highlighter)
      p <- p + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p <- p + geom_point(data=d_order[d_order$SNP %in% highlight_snp, ], aes(x=pos_index, y=pval), colour=highlighter)
    }
  }
  if(!missing(highlight_p)){
    if("Shape" %in% names(d)){
      p <- p + geom_point(data=d_order[d_order$pvalue < highlight_p, ], aes(x=pos_index, y=pval, shape=Shape, frame=Frame), colour=highlighter)
      p <- p + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p <- p + geom_point(data=d_order[d_order$pvalue < highlight_p, ], aes(x=pos_index, y=pval), colour=highlighter)
    }
  }
  #Add title and y axis title
  p <- p + ggtitle(title) + ylab(yaxislab)
  #Add pvalue threshold line
  if(!missing(line)){
    p <- p + geom_hline(yintercept = redline, colour="red")
  }

  #Animate and save
  print(paste("Saving plot to ", file, ".", ext, sep=""))
  ap <- gganimate(p)
  gganimate_save(ap, filename=paste(file, ".", ext, sep=""))
  return(ap)
}
