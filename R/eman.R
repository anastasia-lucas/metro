#' eman
#'
#' Create Manhattan plots for EWAS
#' Note: There is an issue with dev.off() if using RStudio
#' Dependencies: ggplot2
#' Suggested: RColorBrewer
#' @param d data frame, columns one and two must be Variable and pvalue, Shape and Color optional
#' @param format default plotman, other CLARITE and PLATO
#' @param groups optional file containing variable grouping information
#' @param line optional pvalue threshold to draw red line at
#' @param log10 plot -log10() of pvalue column, boolean
#' @param yaxis label for y-axis, automatically set if log10=TRUE
#' @param opacity opacity of points, from 0 to 1, useful for dense plots
#' @param title optional string for plot title
#' @param color1 first alternating color
#' @param color2 second alternating color
#' @param annotate_var list of variables to annotate
#' @param annotate_p pvalue threshold to annotate
#' @param highlight_var list of variables to highlight
#' @param highlight_p pvalue threshold to highlight
#' @param highlighter color to highlight
#' @param groupcolors named list of colors for data in 'Color' column
#' @param file file name of saved image
#' @param hgt height of plot in inches
#' @param wi width of plot in inches
#' @param res resolution of plot in pixels per inch
#' @return png image(s)
#' @export
#' @examples
#' eman(d, format, groups, line, title=NULL, morecolors=FALSE, file="eman", hgt=7, wi=12, res=300 )

eman <- function(d, format="plotman", groups, line, log10=TRUE, yaxis, opacity=1, title=NULL, annotate_var, annotate_p, highlight_var, highlight_p, highlighter="red", color1="#AAAAAA", color2="#4D4D4D", groupcolors, file="eman", hgt=7, wi=12, res=300){
  if (!requireNamespace(c("ggplot2"), quietly = TRUE)==TRUE) {
    stop("Please install ggplot2 to create visualization.", call. = FALSE)
  } else {
    require("ggplot2")
  }

  if(format=="CLARITE"){
  #  if("Variable_pvalue" %in% names(d) & "LRT_pvalue" %in% names(d)){
  #    d$pvalue <- ifelse(!is.na(d$Variable_pvalue), d$Variable_pvalue, ifelse(!is.na(d$LRT_pvalue), d$LRT_pvalue, NA))
  #  } else if("Variable_pvalue" %in% names(d)){
  #    d$pvalue <- d$Variable_pvalue
  #  } else if("LRT_pvalue" %in% names(d)){
  #    d$pvalue <- d$LRT_pvalue
  #  }
    stop("CLARITE format coming soon...")
  } else if(format=="PLATO"){
    stop("PLATO format coming soon...")
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

  if(missing(groups)){
    if("Shape" %in% names(d)){
      p <- ggplot() + geom_point(data=d, aes(x=factor(Variable), y=-log10(pvalue), shape=factor(Shape)), alpha=opacity)
      p <- p + theme(axis.text.x = element_text(angle=90), axis.title.x=element_blank(), legend.position="bottom", legend.title=element_blank())
    } else {
      p <- ggplot(d, aes(x=factor(Variable), y=-log10(pvalue))) + geom_point() + theme(axis.text.x = element_text(angle=90), axis.title.x=element_blank())
    }
  } else {
    subd <- d[, colnames(d) %in% c("Variable", "pvalue", "Shape")]
    colnames(groups)[2] <- "Color"
    dg <- merge(subd, groups, by="Variable")
    dg$Color <- factor(dg$Color)

    #Order variables according to group
    dg_order <- dg[order(dg$Color, dg$Variable), ]
    dg_order$pos_index <- seq.int(nrow(dg_order))

    #Set up dataframe with color and position info
    maxRows <- by(dg_order, dg_order$Color, function(x) x[which.max(x$pos_index),])
    minRows <- by(dg_order, dg_order$Color, function(x) x[which.min(x$pos_index),])
    milimits <- do.call(rbind, minRows)
    malimits <- do.call(rbind, maxRows)
    lims <- merge(milimits, malimits, by="Color")
    if("Shape" %in% names(dg)){
      names(lims) <- c("Color", "Varx", "px", "Shapex", "posmin", "Vary", "py", "Shapey", "posmax")
    } else {
      names(lims) <- c("Color", "Varx", "px", "posmin", "Vary", "py", "posmax")
    }
    lims$av <- (lims$posmin + lims$posmax)/2
    lims$shademap <- rep(c("shade_ebebeb","shade_fffff"), each=1)

    #Set up color palette
    ####FIX
    ncolors <- nlevels(lims$Color)
    if(morecolors==TRUE){
      if (!requireNamespace(c("RColorBrewer"), quietly = TRUE)==TRUE) {
        stop("Please install RColorBrewer to add color attribute.", call. = FALSE)
      }else {
        require("RColorBrewer")
      }
      getPalette = colorRampPalette(brewer.pal(11, "Spectral"))
      newcols <- c(getPalette(ncolors), "#EBEBEB", "#FFFFFF")
    } else {
      newcols <-c(rep(x=c("#53868B", "#4D4D4D"), length.out=ncolors, each=1), "#EBEBEB", "#FFFFFF")
    }
    names(newcols) <-c(levels(factor(lims$Color)), levels(factor(lims$shademap)))

    #Start plotting
    p <- ggplot() + geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = 0, ymax = Inf, fill=factor(shademap)), alpha = 0.5)
    #Add shape info if available
    if("Shape" %in% names(dg)){
      p <- p + geom_point(data=dg_order, aes(x=pos_index, y=-log10(pvalue), color=Color, shape=factor(Shape)), alpha=opacity)
    } else {
      p <- p + geom_point(data=dg_order, aes(x=pos_index, y=-log10(pvalue), color=Color), alpha=opacity)
    }
    p <- p + scale_x_continuous(breaks=lims$av, labels=lims$Color, expand=c(0,0))
    p <- p + geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = -Inf, ymax = 0, fill=Color), alpha = 1)
    p <- p + scale_colour_manual(name = "Color",values = newcols, guides(alpha=FALSE)) + scale_fill_manual(name = "Color",values = newcols, guides(alpha=FALSE))
    p <- p + theme(axis.text.x=element_text(angle=90), panel.grid.minor.x = element_blank(), panel.grid.major.x=element_blank(), axis.title.x=element_blank(), legend.position="bottom", legend.title=element_blank())
  }

  #Add title and y axis title
  p <- p + ggtitle(title) + ylab(expression(paste("-log"[10], "(p-value)", sep="")))

  #Add pvalue threshold line
  if(!missing(line)){
    p <- p + geom_hline(yintercept = -log10(line), colour="red")
  }
  print(paste("Saving plot to ", file, ".png", sep=""))
  ggsave(p, filename=paste(file, ".png", sep=""), dpi=res, units="in", height=hgt, width=wi)
  print(p)

  return(p)
}

