#' eman
#'
#' Create Manhattan plots for EWAS
#' Note: There is an issue with dev.off() if using RStudio
#' Dependencies: ggplot2
#' Suggested: RColorBrewer
#' @param d data frame, columns one and two must be Variable and pvalue, Group, Shape and Color optional
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

eman <- function(d, line, log10=TRUE, yaxis, opacity=1, title=NULL, annotate_var, annotate_p, highlight_var, highlight_p, highlighter="red", color1="#AAAAAA", color2="#4D4D4D", groupcolors, file="eman", hgt=7, wi=12, res=300){
  if (!requireNamespace(c("ggplot2"), quietly = TRUE)==TRUE) {
    stop("Please install ggplot2 to create visualization.", call. = FALSE)
  } else {
    require("ggplot2")
  }

  d$rowid <- seq.int(nrow(d))

  #Info for y-axis
  if(log10==TRUE){
    d$pval <- -log10(d$pvalue)
    yaxislab <- expression(paste("-log"[10], "(p-value)", sep=""))
    if(!missing(line)) {redline <- -log10(line)}
  } else {
    d$pval <- d$pvalue
    yaxislab <- yaxis
    if(!missing(line)) {redline <- line}
  }

  #Add if color here
  if(!"Group" %in% colnames(d)){
    d_order <- merge(d, dpval, by="rowid")
    dpval <- d_order[, colnames(d_order) %in% c("rowid", "pval"), drop=FALSE]
    if("Shape" %in% names(d)){
      if("Color" %in% names(d)){
        p <- ggplot() + geom_point(data=d_order, aes(x=factor(Variable), y=pval, shape=factor(Shape), color=Color), alpha=opacity)
      } else {
        p <- ggplot() + geom_point(data=d_order, aes(x=factor(Variable), y=pval, shape=factor(Shape)), alpha=opacity)
      }
      p <- p + theme(axis.text.x = element_text(angle=90), axis.title.x=element_blank(), legend.position="bottom", legend.title=element_blank())
    } else {
      if("Color" %in% names(d)){
        p <- ggplot(d_order, aes(x=factor(Variable), y=pval, color=Color)) + geom_point()
        p <- p + theme(axis.text.x = element_text(angle=90), axis.title.x=element_blank(), legend.position="bottom", legend.title=element_blank())
      } else {
        p <- ggplot(d_order, aes(x=factor(Variable), y=pval)) + geom_point() + theme(axis.text.x = element_text(angle=90), axis.title.x=element_blank())
      }
    }
  } else {
    #Create position index
    subd <- d[, colnames(d) %in% c("Variable", "Group", "pvalue", "rowid")]
    d_order <- subd[order(subd$Group, subd$Variable),]
    d_order$pos_index <- seq.int(nrow(d_order))


    #Set up dataframe with position info
    subd <- d_order[, colnames(d_order)!="rowid"]
    maxRows <- by(subd, subd$Group, function(x) x[which.max(x$pos_index),])
    minRows <- by(subd, subd$Group, function(x) x[which.min(x$pos_index),])
    milimits <- do.call(rbind, minRows)
    malimits <- do.call(rbind, maxRows)
    lims <- merge(milimits, malimits, by="Group")
    names(lims) <- c("Color", "Varx", "px", "posmin", "Vary", "py", "posmax")
    lims$av <- (lims$posmin + lims$posmax)/2
    lims$shademap <- rep(c("shade_ffffff","shade_ebebeb"), each=1)

    #Set up colors
    nvarcolors <- nlevels(factor(lims$Color))
    if("Color" %in% names(d)){
      if(!missing(groupcolors)){
        dcols <- c(rep(x=c(color1, color2), length.out=nvarcolors, each=1), "#FFFFFF", "#EBEBEB")
        names(dcols) <-c(levels(factor(lims$Color)), "shade_ffffff", "shade_ebebeb")
        newcols <- c(dcols, groupcolors)
      } else {
        if (!requireNamespace(c("RColorBrewer"), quietly = TRUE)==TRUE) {
          stop("Please install RColorBrewer to add color attribute.", call. = FALSE)
        } else {
          require("RColorBrewer", quietly=TRUE)
        }
        ngroupcolors <- nlevels(factor(d$Color))
        getPalette = colorRampPalette(brewer.pal(11, "Spectral"))
        newcols <- c(rep(x=c(color1, color2), length.out=nvarcolors, each=1), getPalette(ngroupcolors), "#FFFFFF", "#EBEBEB")
        names(newcols) <-c(levels(factor(lims$Color)), levels(factor(d$Color)), "shade_ffffff", "shade_ebebeb")
      }
    } else {
      #Color by Group instead
      names(d_order)[names(d_order)=="Group"] <- "Color"
      newcols <-c(rep(x=c(color1, color2), length.out=nvarcolors, each=1), "#FFFFFF", "#EBEBEB")
      names(newcols) <-c(levels(factor(lims$Color)), "shade_ffffff", "shade_ebebeb")
    }

    #Start plotting
    dinfo <- d[, colnames(d) %in% c("rowid", "Color", "Shape", "pval"), drop=FALSE]
    d_order <- merge(d_order, dinfo, by="rowid")
    p <- ggplot() + geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = 0, ymax = Inf, fill=factor(shademap)), alpha = 0.5)
    #Add shape info if available
    if("Shape" %in% names(d)){
      p <- p + geom_point(data=d_order, aes(x=pos_index, y=pval, color=Color, shape=factor(Shape)), alpha=opacity)
    } else {
      p <- p + geom_point(data=d_order, aes(x=pos_index, y=pval, color=Color), alpha=opacity)
    }

    p <- p + scale_x_continuous(breaks=lims$av, labels=lims$Color, expand=c(0,0))
    p <- p + geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = -Inf, ymax = 0, fill=Color), alpha = 1)
    #p <- p + scale_colour_manual(name = "Color",values = newcols, guides(alpha=FALSE)) + scale_fill_manual(name = "Color",values = newcols, guides(alpha=FALSE))
    p <- p + theme(axis.text.x=element_text(angle=90), panel.grid.minor.x = element_blank(), panel.grid.major.x=element_blank(), axis.title.x=element_blank(), legend.position="bottom", legend.title=element_blank())
  }
  if("Color" %in% names(d)){
    #Add legend
    p <- p + scale_colour_manual(name = "Color", values = newcols) + scale_fill_manual(name = "Color", values = newcols, guides(alpha=FALSE))
  } else {
    #Don't
    p <- p + scale_colour_manual(name = "Color", values = newcols, guides(alpha=FALSE)) + scale_fill_manual(name = "Color", values = newcols, guides(alpha=FALSE))
  }
  if(!missing(annotate_p)){
    if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE) {
      print("Consider installing 'ggrepel' for improved text annotation")
      p <- p + geom_text(data=d_order[d_order$pvalue < annotate_p,], aes(pos_index,pval,label=Variable))
    } else {
      require("ggrepel", quietly = TRUE)
      p <- p + geom_text_repel(data=d_order[d_order$pvalue < annotate_p,], aes(pos_index,pval,label=Variable))
    }
  }
  if(!missing(annotate_var)){
    if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE){
      print("Consider installing 'ggrepel' for improved text annotation")
      p <- p + geom_text(data=d_order[d_order$Variable %in% annotate_var,], aes(pos_index,pval,label=Variable))
    } else {
      require("ggrepel", quietly = TRUE)
      p <- p + geom_text_repel(data=d_order[d_order$Variable %in% annotate_var,], aes(pos_index,pval,label=Variable))
    }
  }
  #Highlight if given
  if(!missing(highlight_var)){
    if("Shape" %in% names(d)){
      p <- p + geom_point(data=d_order[d_order$Variable %in% highlight_var, ], aes(x=pos_index, y=pval, shape=Shape), colour=highlighter)
      p <- p + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p <- p + geom_point(data=d_order[d_order$Variable %in% highlight_var, ], aes(x=pos_index, y=pval), colour=highlighter)
    }
  }
  if(!missing(highlight_p)){
    if("Shape" %in% names(d)){
      p <- p + geom_point(data=d_order[d_order$pvalue < highlight_p, ], aes(x=pos_index, y=pval, shape=Shape), colour=highlighter)
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

  #Save
  print(paste("Saving plot to ", file, ".png", sep=""))
  ggsave(p, filename=paste(file, ".png", sep=""), dpi=res, units="in", height=hgt, width=wi)
  print(p)

  return(p)
}

