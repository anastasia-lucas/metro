#' igman
#'
#' Create Interactive Manhattan plots for GWAS
#' Dependencies: ggplot2, ggiraph
#' Suggested: RColorBrewer
#' @param d data frame, if not plato or plink format, must contain SNP, CHR, POS, pvalue columns, optional Shape, Color, and Info
#' @param format format of input
#' @param line optional pvalue threshold to draw red line at
#' @param log10 plot -log10() of pvalue column, boolean
#' @param yaxis label for y-axis, automatically set if log10=TRUE
#' @param title optional string for plot title
#' @param chrcolor1 first alternating color for chromosome
#' @param chrcolor2 second alternating color for chromosome
#' @param highlight_snp list of snps to highlight
#' @param highlight_p pvalue threshold to highlight
#' @param highlighter color to highlight
#' @param groupcolors named list of colors for data in 'Color' column
#' @param db choose database to connect to (dbSNP or GWAScatalog)
#' @param moreinfo includes more information on hover, refers to Info column
#' @param file file name of saved image
#' @param hgt height of plot in inches
#' @param wi width of plot in inches
#' @return html file
#' @export
#' @examples
#' igman(d, format, line, log10, yaxis, title, chrcolor1, chrcolor2, groupcolors, db, moreinfo, file, hgt, wi)

igman <- function(d, format="plotman", line, log10=TRUE, yaxis, highlight_snp, highlight_p, highlighter="red", title=NULL, chrcolor1="#AAAAAA", chrcolor2="#4D4D4D", groupcolors, db, moreinfo=FALSE, file="igman", hgt=7, wi=12){
  if (!requireNamespace(c("ggplot2"), quietly = TRUE)==TRUE|!requireNamespace(c("ggiraph"), quietly = TRUE)==TRUE) {
    stop("Please install ggplot2 and ggiraph to create interactive visualization.", call. = FALSE)
  } else {
    require("ggplot2", quietly=TRUE)
    require("ggiraph", quietly=TRUE)
  }

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
  d_order <- d[order(d$CHR, d$POS), ]
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
  lims$shademap <- rep(c("shade_ffffff","shade_ebebeb"), each=1)


  #Set up tooltip
  ###See what info would be useful here i.e. SNP or something else
  d_order$tooltip <- if (moreinfo==TRUE) c(paste0(d_order$SNP, "\n Add'l:", d_order$Info)) else d_order$SNP

  #Set up onclick
  if(!missing(db)){
    if(db=="GWAScatalog"){
      d_order$onclick <- sprintf("window.open(\"%s%s\")","https://www.ebi.ac.uk/gwas/search?query=", as.character(d_order$SNP))
    } else if(db=="dbSNP"){
      d_order$onclick <- sprintf("window.open(\"%s%s\")","https://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?searchType=adhoc_search&type=rs&rs=", as.character(d_order$SNP))
    }
  } else {
    d_order$onclick <- NA
  }

  #Set up colors
  nchrcolors <- nlevels(factor(lims$Color))
  if("Color" %in% names(d)){
    if(!missing(groupcolors)){
      dcols <- c(rep(x=c(chrcolor1, chrcolor2), length.out=nchrcolors, each=1), "#FFFFFF", "#EBEBEB")
      names(dcols) <-c(levels(factor(lims$Color)), "shade_fffff", "shade_ebebeb")
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
      names(newcols) <-c(levels(factor(lims$Color)), levels(factor(d_order$Color)), "shade_fffff", "shade_ebebeb")
    }
  } else {
    #Color by CHR instead
    colnames(d_order)[2] <- "Color"
    newcols <-c(rep(x=c(chrcolor1, chrcolor2), length.out=nchrcolors, each=1), "#FFFFFF", "#EBEBEB")
    names(newcols) <-c(levels(factor(lims$Color)), "shade_ffffff", "shade_ebebeb")
  }

  #Info for y-axis
  if(log10==TRUE){
    d_order$pval <- -log10(d_order$pvalue)
    yaxislab <- expression(paste("-log"[10], "(p-value)", sep=""))
    redline <- -log10(line)
  } else {
    d_order$pval <- d_order$pvalue
    yaxislab <- yaxis
    redline <- line
  }

  #Start plotting
  p <- ggplot() + geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = 0, ymax = Inf, fill=factor(shademap)), alpha = 0.5)
  #Add shape info if available
  if("Shape" %in% names(d)){
    p <- p + geom_point_interactive(data=d_order, aes(x=pos_index, y=pval, tooltip=tooltip, onclick=onclick, color=factor(Color), shape=factor(Shape)))
  } else {
    p <- p + geom_point_interactive(data=d_order, aes(x=pos_index, y=pval, tooltip=tooltip, onclick=onclick, color=factor(Color)))
  }
  p <- p + scale_x_continuous(breaks=lims$av, labels=lims$Color, expand=c(0,0))
  p <- p + geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = -Inf, ymax = 0, fill=as.factor(Color)), alpha = 1)
  if("Color" %in% names(d)){
    #Add legend
    p <- p + scale_colour_manual(name = "Color", values = newcols) + scale_fill_manual(name = "Color",values = newcols, guides(alpha=FALSE))
  } else {
    #Don't
    p <- p + scale_colour_manual(name = "Color", values = newcols, guides(alpha=FALSE)) + scale_fill_manual(name = "Color", values = newcols, guides(alpha=FALSE))
  }
  p <- p + theme(axis.text.x=element_text(angle=90), panel.grid.minor.x = element_blank(), panel.grid.major.x=element_blank(), axis.title.x=element_blank(), legend.position="bottom", legend.title=element_blank())
  #Highlight if given
  if(!missing(highlight_snp)){
    if("Shape" %in% names(d)){
      p <- p + geom_point_interactive(data=d_order[d_order$SNP %in% highlight_snp, ], aes(x=pos_index, y=pval, shape=Shape, tooltip=tooltip, onclick=onclick), colour=highlighter)
      p <- p + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p <- p + geom_point_interactive(data=d_order[d_order$SNP %in% highlight_snp, ], aes(x=pos_index, y=pval, tooltip=SNP, onclick=onclick), colour=highlighter)
    }
  }
  if(!missing(highlight_p)){
    if("Shape" %in% names(d)){
      p <- p + geom_point_interactive(data=d_order[d_order$pvalue < highlight_p, ], aes(x=pos_index, y=pval, shape=Shape, tooltip=tooltip, onclick=onclick), colour=highlighter)
      p <- p + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p <- p + geom_point_interactive(data=d_order[d_order$pvalue < highlight_p, ], aes(x=pos_index, y=pval, tooltip=tooltip, onclick=onclick), colour=highlighter)
    }
  }
  #Add title and y axis title
  p <- p + ggtitle(title) + ylab(expression(paste("-log"[10], "(p-value)", sep="")))
  #Add pvalue threshold line
  if(!missing(line)){
    p <- p + geom_hline(yintercept = redline, colour="red")
  }

  #Save
  print(paste("Saving plot to ", file, ".html", sep=""))
  tooltip_css <- "background-color:black;color:white;padding:6px;border-radius:15px 15px 15px 15px;"
  ip <- ggiraph(code=print(p), tooltip_extra_css = tooltip_css, tooltip_opacity = 0.75, zoom_max = 6)
  htmlwidgets::saveWidget(widget=ip, file=paste(file, ".html"))
  return(ip)

}
