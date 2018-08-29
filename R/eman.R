#' eman
#'
#' Create Manhattan plots for EWAS
#' @param d data frame, columns are Variable, pvalue, and Group; Shape and Color optional
#' @param line optional pvalue threshold to draw red line at
#' @param log10 plot -log10() of pvalue column, boolean
#' @param yaxis label for y-axis, automatically set if log10=TRUE
#' @param opacity opacity of points, from 0 to 1, useful for dense plots
#' @param title optional string for plot title
#' @param color1 first alternating color
#' @param color2 second alternating color
#' @param annotate_var vector of variables to annotate
#' @param annotate_p pvalue threshold to annotate
#' @param highlight_var vector of variables to highlight
#' @param highlight_p pvalue threshold to highlight
#' @param highlighter color to highlight
#' @param groupcolors named vector of colors for data in 'Color' column
#' @param background variegated or white
#' @param grpblocks boolean, turns on x-axis group marker blocks
#' @param file file name of saved image
#' @param hgt height of plot in inches
#' @param wi width of plot in inches
#' @param res resolution of plot in pixels per inch
#' @import ggplot2
#' @return png image(s)
#' @export
#' @family EWAS functions
#' @family static plotting functions
#' @seealso \code{\link{ieman}}, \code{\link{aeman}}, \code{\link{gman}}, \code{\link{pheman}}
#' @examples
#' data(ewas)
#' eman(d=ewas, title="EWAS", line=0.001, annotate_p=0.001, color1="#A23B72", color2="#2A84AA",
#' highlight_p=0.001, highlighter="green")

eman <- function(d, line, log10=TRUE, yaxis, opacity=1, title=NULL, annotate_var, annotate_p, highlight_var, highlight_p, highlighter="red", color1="#AAAAAA", color2="#4D4D4D", groupcolors, background="variegated", grpblocks=FALSE, file="eman", hgt=7, wi=12, res=300){
  if (!requireNamespace(c("ggplot2"), quietly = TRUE)==TRUE) {
    stop("Please install ggplot2 to create visualization.", call. = FALSE)
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

  #Theme options
  yaxismin <- min(d$pval)
  backpanel <- ifelse(background=="white", "NULL", "geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = yaxismin, ymax = Inf, fill=factor(shademap)), alpha = 0.5)" )

  #Allow more than 6 shapes
  if("Shape" %in% names(d)){
    allshapes <- c(16,15,17,3,7,8,0:2,4:6,9:14,18:25,33:127)
    shapevector <- allshapes[1:nlevels(as.factor(d$Shape))]
  }

  #Save to merge later
  d$rowid <- seq.int(nrow(d))
  dinfo <- d[, colnames(d) %in% c("rowid", "Color", "Shape", "pval"), drop=FALSE]

  #If no group, plot raw data
  #if(!"Group" %in% colnames(d)){
  #  d_order <- merge(d, dinfo, by="rowid")
  #  if("Shape" %in% names(d)){
  #    if("Color" %in% names(d)){
  #      p <- ggplot() + geom_point(data=d_order, aes(x=factor(Variable), y=pval, shape=factor(Shape), color=Color), alpha=opacity) + scale_shape_manual(values=shapevector)
  #    } else {
  #      p <- ggplot() + geom_point(data=d_order, aes(x=factor(Variable), y=pval, shape=factor(Shape)), alpha=opacity) + scale_shape_manual(values=shapevector)
  #    }
  #    p <- p + theme(axis.text.x = element_text(angle=90), axis.title.x=element_blank(), legend.position="bottom", legend.title=element_blank())
  #  } else {
  #    if("Color" %in% names(d)){
  #      p <- ggplot(d_order, aes(x=factor(Variable), y=pval, color=Color)) + geom_point()
  #      p <- p + theme(axis.text.x = element_text(angle=90), axis.title.x=element_blank(), legend.position="bottom", legend.title=element_blank())
  #    } else {
  #      p <- ggplot(d_order, aes(x=factor(Variable), y=pval)) + geom_point() + theme(axis.text.x = element_text(angle=90), axis.title.x=element_blank())
  #    }
  #  }
  #} else {
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
        if(ngroupcolors > 15){
          getPalette = colorRampPalette(brewer.pal(11, "Spectral"))
          newcols <- c(rep(x=c(color1, color2), length.out=nvarcolors, each=1), getPalette(ngroupcolors), "#FFFFFF", "#EBEBEB")
        } else {
          pal <- pal <- c("#009292", "#920000", "#490092", "#db6d00", "#24ff24",
                          "#ffff6d", "#000000", "#006ddb", "#004949","#924900",
                          "#ff6db6", "#6db6ff","#b66dff", "#ffb6db","#b6dbff")
          newcols <- c(rep(x=c(color1, color2), length.out=nvarcolors, each=1),pal[1:ngroupcolors], "#FFFFFF", "#EBEBEB")
        }
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
    p <- ggplot() + eval(parse(text=backpanel))
    #Add shape info if available
    if("Shape" %in% names(d)){
      p <- p + geom_point(data=d_order, aes(x=pos_index, y=pval, color=Color, shape=factor(Shape)), alpha=opacity) + scale_shape_manual(values=shapevector)
    } else {
      p <- p + geom_point(data=d_order, aes(x=pos_index, y=pval, color=Color), alpha=opacity)
    }
    p <- p + scale_x_continuous(breaks=lims$av, labels=lims$Color, expand=c(0,0))
    if(grpblocks==TRUE){p <- p + geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = -Inf, ymax = min(d_order$pval), fill=as.factor(Color)), alpha = 1)}
    #p <- p + scale_colour_manual(name = "Color",values = newcols, guides(alpha=FALSE)) + scale_fill_manual(name = "Color",values = newcols, guides(alpha=FALSE))
    p <- p + theme(axis.text.x=element_text(angle=90), panel.grid.minor.x = element_blank(), panel.grid.major.x=element_blank(), axis.title.x=element_blank(), legend.position="bottom", legend.title=element_blank())
  #}
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
      p <- p + ggrepel::geom_text_repel(data=d_order[d_order$pvalue < annotate_p,], aes(pos_index,pval,label=Variable))
    }
  }
  if(!missing(annotate_var)){
    if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE){
      print("Consider installing 'ggrepel' for improved text annotation")
      p <- p + geom_text(data=d_order[d_order$Variable %in% annotate_var,], aes(pos_index,pval,label=Variable))
    } else {
      p <- p + ggrepel::geom_text_repel(data=d_order[d_order$Variable %in% annotate_var,], aes(pos_index,pval,label=Variable))
    }
  }
  #Highlight if given
  if(!missing(highlight_var)){
    if("Shape" %in% names(d)){
      p <- p + geom_point(data=d_order[d_order$Variable %in% highlight_var, ], aes(x=pos_index, y=pval, shape=Shape), colour=highlighter) + scale_shape_manual(values=shapevector)
      p <- p + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p <- p + geom_point(data=d_order[d_order$Variable %in% highlight_var, ], aes(x=pos_index, y=pval), colour=highlighter)
    }
  }
  if(!missing(highlight_p)){
    if("Shape" %in% names(d)){
      p <- p + geom_point(data=d_order[d_order$pvalue < highlight_p, ], aes(x=pos_index, y=pval, shape=Shape), colour=highlighter) + scale_shape_manual(values=shapevector)
      p <- p + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p <- p + geom_point(data=d_order[d_order$pvalue < highlight_p, ], aes(x=pos_index, y=pval), colour=highlighter)
    }
  }

  #Add title and y axis title
  p <- p + ggtitle(title) + ylab(yaxislab)

  #Add pvalue threshold line
  if(!missing(line)){p <- p + geom_hline(yintercept = redline, colour="red")}

  if(grpblocks==TRUE){
    p <- p+ylim(c(yaxismin,max(d_order$pval)))
  } else {
    p <- p+scale_y_continuous(limits=c(yaxismin, max(d_order$pval)),expand=expand_scale(mult=c(0,0.1)))
  }
  if(background=="white"){p <- p + theme(panel.background = element_rect(fill="white"))}

  #Save
  print(paste("Saving plot to ", file, ".png", sep=""))
  ggsave(p, filename=paste(file, ".png", sep=""), dpi=res, units="in", height=hgt, width=wi)
  print(p)

  return(p)
}

