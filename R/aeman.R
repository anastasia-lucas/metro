#' aeman
#'
#' Create animated Manhattan plots for EWAS
#' Note: There is an issue with dev.off() if using RStudio
#' Dependencies: ggplot2
#' Suggested: RColorBrewer
#' @param d data frame, columns one, two, and three must be Variable, pvalue, and Frame; Group, Shape and Color optional
#' @param line optional pvalue threshold to draw red line at
#' @param log10 plot -log10() of pvalue column, boolean
#' @param yaxis label for y-axis, automatically set if log10=TRUE
#' @param opacity opacity of points, from 0 to 1, useful for dense plots
#' @param title optional string for plot title
#' @param color1 first alternating color
#' @param color2 second alternating color
#' @param highlight_var list of variables to highlight
#' @param highlight_p pvalue threshold to highlight
#' @param highlighter color to highlight
#' @param groupcolors named list of colors for data in 'Color' column
#' @param file file name of saved image
#' @param ext file type to save, "gif" or "html"
#' @param hgt height of plot in pixels
#' @param wi width of plot in pixels
#' @return .gif or .html file
#' @export
#' @examples
#' aeman(d, format, groups, line, title=NULL, file="eman", hgt=1300, wi=800, )

aeman <- function(d, line, log10=TRUE, yaxis, opacity=1, title=NULL, highlight_var, highlight_p, highlighter="red", color1="#AAAAAA", color2="#4D4D4D", groupcolors, file="aeman", hgt=800, wi=1300){
  if (!requireNamespace(c("ggplot2"), quietly = TRUE)==TRUE) {
    stop("Please install ggplot2 to create visualization.", call. = FALSE)
  } else {
    require("ggplot2")
  }

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

  #Save to merge later
  d$rowid <- seq.int(nrow(d))
  dinfo <- d[, colnames(d) %in% c("rowid", "Color", "Shape", "pval"), drop=FALSE]

  #If no group, plot raw data
  if(!"Group" %in% colnames(d)){
    d_order <- merge(d, dinfo, by="rowid")
    if("Shape" %in% names(d)){
      if("Color" %in% names(d)){
        p <- ggplot() + geom_point(data=d_order, aes(x=factor(Variable), y=pval, shape=Shape, color=Color, frame=Frame), alpha=opacity)
      } else {
        p <- ggplot() + geom_point(data=d_order, aes(x=factor(Variable), y=pval, shape=Shape, frame=Frame), alpha=opacity)
      }
      p <- p + theme(axis.text.x = element_text(angle=90), axis.title.x=element_blank(), legend.position="bottom", legend.title=element_blank())
    } else {
      if("Color" %in% names(d)){
        p <- ggplot(d_order, aes(x=factor(Variable), y=pval, color=Color, frame=Frame)) + geom_point()
        p <- p + theme(axis.text.x = element_text(angle=90), axis.title.x=element_blank(), legend.position="bottom", legend.title=element_blank())
      } else {
        p <- ggplot(d_order, aes(x=factor(Variable), y=pval, frame=Frame)) + geom_point() + theme(axis.text.x = element_text(angle=90), axis.title.x=element_blank())
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
    lims$shademap <- rep(c("shade_ffffff","shade_ebebeb"), each=1, length.out=nrow(lims))

    #Set up colors
    nvarcolors <- nlevels(factor(lims$Color))
    if("Color" %in% names(d)){
      #Color by Color column
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
    d_order <- merge(d_order, dinfo, by="rowid")
    p <- ggplot() + geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = 0, ymax = Inf, fill=factor(shademap)), alpha = 0.5)
    #Add shape info if available
    if("Shape" %in% names(d)){
      p <- p + geom_point(data=d_order, aes(x=pos_index, y=pval, color=Color, shape=Shape, frame=Frame), alpha=opacity)
    } else {
      p <- p + geom_point(data=d_order, aes(x=pos_index, y=pval, color=Color, frame=Frame), alpha=opacity)
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
  #if(!missing(annotate_p)){
  #  if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE) {
  #    print("Consider installing 'ggrepel' for improved text annotation")
  #    p <- p + geom_text(data=d_order[d_order$pvalue < annotate_p,], aes(pos_index,pval,label=Variable))
  #  } else {
  #    require("ggrepel", quietly = TRUE)
  #    p <- p + geom_text_repel(data=d_order[d_order$pvalue < annotate_p,], aes(pos_index,pval,label=Variable))
  #  }
  #}
  #if(!missing(annotate_var)){
  #  if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE){
  #    print("Consider installing 'ggrepel' for improved text annotation")
  #    p <- p + geom_text(data=d_order[d_order$Variable %in% annotate_var,], aes(pos_index,pval,label=Variable))
  #  } else {
  #    require("ggrepel", quietly = TRUE)
  #    p <- p + geom_text_repel(data=d_order[d_order$Variable %in% annotate_var,], aes(pos_index,pval,label=Variable))
  #  }
  #}
  #Highlight if given
  if(!missing(highlight_var)){
    if("Shape" %in% names(d)){
      p <- p + geom_point(data=d_order[d_order$Variable %in% highlight_var, ], aes(x=pos_index, y=pval, shape=Shape, frame=Frame), colour=highlighter)
      p <- p + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p <- p + geom_point(data=d_order[d_order$Variable %in% highlight_var, ], aes(x=pos_index, y=pval, frame=Frame), colour=highlighter)
    }
  }
  if(!missing(highlight_p)){
    if("Shape" %in% names(d)){
      p <- p + geom_point(data=d_order[d_order$pvalue < highlight_p, ], aes(x=pos_index, y=pval, shape=Shape, frame=Frame), colour=highlighter)
      p <- p + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p <- p + geom_point(data=d_order[d_order$pvalue < highlight_p, ], aes(x=pos_index, y=pval, frame=Frame), colour=highlighter)
    }
  }

  #Add title and y axis title
  p <- p + ggtitle(title) + ylab(yaxislab)

  #Add pvalue threshold line
  if(!missing(line)){
    p <- p + geom_hline(yintercept = redline, colour="red")
  }

  #Save
  print(paste("Saving plot to ", file, ".", ext, sep=""))
  ap <- gganimate(p)
  gganimate_save(ap, filename=paste(file, ".", ext, sep=""), ani.height=hgt, ani.width=wi)
  return(ap)
}

