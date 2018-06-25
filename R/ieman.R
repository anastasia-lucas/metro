#' ieman
#'
#' Create interactive Manhattan plots for EWAS
#' Note: There is an issue with dev.off() if using RStudio
#' Dependencies: ggplot2
#' Suggested: RColorBrewer
#' @param d data frame, columns one and two must be Variable and pvalue; Group, Shape, Color, and Info optional
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
#' @param db query address, ex. "https://www.google.com/search?q"
#' @param moreinfo includes more information on hover, refers to Info column
#' @param file file name of saved image
#' @param hgt height of plot in inches
#' @param wi width of plot in inches
#' @param res resolution of plot in pixels per inch
#' @return png image(s)
#' @export
#' @examples
#' ieman(d, groups, line, title=NULL, morecolors=FALSE, file="eman", hgt=7, wi=12, res=300 )

ieman <- function(d, line, log10=TRUE, yaxis, opacity=1, title=NULL, highlight_var, highlight_p, highlighter="red", color1="#AAAAAA", color2="#4D4D4D", groupcolors, db, moreinfo=FALSE, file="ieman", hgt=7, wi=12, res=300){
  if (!requireNamespace(c("ggplot2"), quietly = TRUE)==TRUE|!requireNamespace(c("ggiraph"), quietly = TRUE)==TRUE) {
    stop("Please install ggplot2 and ggiraph to create interactive visualization.", call. = FALSE)
  } else {
    require("ggplot2", quietly=TRUE)
    require("ggiraph", quietly=TRUE)
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

  #Allow more than 6 shapes
  #3, 4 and 7 to 14 are composite symbols- incompatible with ggiraph
  if("Shape" %in% names(d)){
    allshapes <- c(16,15,17,18,0:2,5:6,19:25,33:127)
    shapevector <- allshapes[1:nlevels(as.factor(d$Shape))]
  }

  #Set up tooltip
  d$tooltip <- if (moreinfo==TRUE) c(paste0(d$Variable, "\n ", d$Info)) else d$Variable

  #Set up onclick
  if(!missing(db)){
    d$onclick <- sprintf("window.open(\"%s%s\")", db, as.character(d$Variable))
  } else {
    d$onclick <- NA
  }

  #Save to merge later
  d$rowid <- seq.int(nrow(d))
  dinfo <- d[, colnames(d) %in% c("rowid", "Color", "Shape", "pval", "tooltip", "onclick"), drop=FALSE]

  #If no group, plot raw data
  if(!"Group" %in% colnames(d)){
    d_order <- merge(d, dinfo, by="rowid")
    if("Shape" %in% names(d)){
      if("Color" %in% names(d)){
        p <- ggplot() + geom_point_interactive(data=d_order, aes(x=factor(Variable), y=pval, shape=factor(Shape), color=Color, onclick=onclick, tooltip=tooltip), alpha=opacity) + scale_shape_manual(values=shapevector)
      } else {
        p <- ggplot() + geom_point_interactive(data=d_order, aes(x=factor(Variable), y=pval, shape=factor(Shape), onclick=onclick, tooltip=tooltip), alpha=opacity) + scale_shape_manual(values=shapevector)
      }
      p <- p + theme(axis.text.x = element_text(angle=90), axis.title.x=element_blank(), legend.position="bottom", legend.title=element_blank())
    } else {
      if("Color" %in% names(d)){
        p <- ggplot() + geom_point_interactive(data=d_order, aes(x=factor(Variable), y=pval, color=Color, tooltip=tooltip, onclick=onclick))
        p <- p + theme(axis.text.x = element_text(angle=90), axis.title.x=element_blank(), legend.position="bottom", legend.title=element_blank())
      } else {
        p <- ggplot() + geom_point_interactive(data=d_order, aes(x=factor(Variable), y=pval, tooltip=tooltip, onclick=onclick)) + theme(axis.text.x = element_text(angle=90), axis.title.x=element_blank())
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
      p <- p + geom_point_interactive(data=d_order, aes(x=pos_index, y=pval, color=Color, shape=factor(Shape), onclick=onclick, tooltip=tooltip), alpha=opacity) + scale_shape_manual(values=shapevector)
    } else {
      p <- p + geom_point_interactive(data=d_order, aes(x=pos_index, y=pval, color=Color, onclick=onclick, tooltip=tooltip), alpha=opacity)
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

  #Highlight if given
  if(!missing(highlight_var)){
    if("Shape" %in% names(d)){
      p <- p + geom_point_interactive(data=d_order[d_order$Variable %in% highlight_var, ], aes(x=pos_index, y=pval, shape=Shape, onclick=onclick, tooltip=tooltip), colour=highlighter) + scale_shape_manual(values=shapevector)
      p <- p + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p <- p + geom_point_interactive(data=d_order[d_order$Variable %in% highlight_var, ], aes(x=pos_index, y=pval, onclick=onclick, tooltip=tooltip), colour=highlighter)
    }
  }
  if(!missing(highlight_p)){
    if("Shape" %in% names(d)){
      p <- p + geom_point_interactive(data=d_order[d_order$pvalue < highlight_p, ], aes(x=pos_index, y=pval, shape=Shape, onclick=onclick, tooltip=tooltip), colour=highlighter) + scale_shape_manual(values=shapevector)
      p <- p + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p <- p + geom_point_interactive(data=d_order[d_order$pvalue < highlight_p, ], aes(x=pos_index, y=pval, onclick=onclick, tooltip=tooltip), colour=highlighter)
    }
  }

  #Add title and y axis title
  p <- p + ggtitle(title) + ylab(yaxislab)

  #Add pvalue threshold line
  if(!missing(line)){
    p <- p + geom_hline(yintercept = redline, colour="red")
  }

  #Save
  print(paste("Saving plot to ", file, ".html", sep=""))
  tooltip_css <- "background-color:black;color:white;padding:6px;border-radius:15px 15px 15px 15px;"
  ip <- ggiraph(code=print(p), tooltip_extra_css = tooltip_css, tooltip_opacity = 0.75, zoom_max = 6, width_svg=wi, height_svg=hgt)
  htmlwidgets::saveWidget(widget=ip, file=paste(file, ".html", sep=""))
  return(ip)
}

