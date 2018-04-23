#' agxgman
#'
#' Create animated heatmap plots for SNPxSNP or GxG interaction
#' Dependencies: ggplot2
#' @param d data frame, must contain SNP1, CHR1, POS1, SNP2, CHR2, POS2, pvalue, Frame columns
#' @param me optional main effect data frame containing SNP, CHR, POS, pvalue, Frame columns, optional Shape
#' @param symmetric boolean, should the plot be symmetric
#' @param highlight_p threshold to diverge gradient color scale, default 0.05, set to "off" for no threshold
#' @param legend title for color legend, default "p-value"
#' @param title optional string for plot title
#' @param high color for high values
#' @param low color for low values
#' @param highlight_high if highlight_p given, color for max of highlight range
#' @param highlight_low if highlight_p given, color for min of highlight range
#' @param file file name of saved image
#' @param ext file extension "html" or "gif"
#' @param hgt height of plot in inches
#' @param wi width of plot in inches
#' @return html or gif file
#' @export
#' @examples
#' agxgman(d, me=NULL, symmetric, highlight_p, legend, title, high, low, highlight_high, highlight_low, file, hgt, wi)

agxgman <- function(d, me, symmetric=TRUE, highlight_p=0.05, legend="p-value", title=NULL, high="#02021e", low="blue", highlight_high="yellow", highlight_low="#fffcd3", file="agxgman", ext="html", hgt=7, wi=7.5, res=300){
  if (!requireNamespace(c("ggplot2"), quietly = TRUE)==TRUE|!requireNamespace(c("gganimate"), quietly = TRUE)==TRUE) {
    stop("Please install ggplot2 and ggiraph to create interactive visualization.", call. = FALSE)
  } else {
    require("ggplot2", quietly=TRUE)
    require("gganimate", quietly=TRUE)
  }

  d$CHR1 <- factor(d$CHR1, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"))
  d$CHR2 <- factor(d$CHR2, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"))

  d$Frame <- factor(d$Frame)

   #Initialize vector for gradient
  if(highlight_p!="off"){
    color_vector <- c(highlight_low, highlight_high, low, high)
    #Values should be min, max < highlight value, highlight value, and max
    value_vector <- c(floor(min(d$pvalue, me$pvalue)), max(d$pval[d$pvalue < highlight_p]), highlight_p,ceiling(max(d$pvalue, me$pvalue)))
    breaks_vector <- c(floor(min(d$pvalue, me$pvalue)), max(d$pval[d$pvalue < highlight_p]), highlight_p, ceiling(max(d$pvalue, me$pvalue)))
    label_vector <- c(floor(min(d$pvalue, me$pvalue)), paste("<",highlight_p,sep=""), highlight_p, ceiling(max(d$pvalue, me$pvalue)))
  } else {
    color_vector <- c(low, high)
    #Values should be min, max < highlight value, highlight value, and max
    value_vector <- c(floor(min(d$pvalue, me$pvalue)), ceiling(max(d$pvalue, me$pvalue)))
    breaks_vector <- c(floor(min(d$pvalue, me$pvalue)), ceiling(max(d$pvalue, me$pvalue)))
    label_vector <- c(floor(min(d$pvalue, me$pvalue)), ceiling(max(d$pvalue, me$pvalue)))
  }

  #Non-symmetric
  if(symmetric==FALSE){
    #Order and assing position by CHR, POS
    d_order <- d[order(d$CHR1, d$POS1), ]
    snp1pos <- unique(d_order[, colnames(d_order) %in% c("SNP1", "POS1", "CHR1")])
    snp1pos$pos_index1 <- seq.int(nrow(snp1pos))
    snp1merge <- snp1pos[, colnames(snp1pos) %in% c("SNP1", "pos_index1")]
    d_order1 <- merge(d_order, snp1merge, by="SNP1")
    d_order <- d[order(d$CHR2, d$POS2), ]
    snp2pos <- unique(d_order[, colnames(d_order) %in% c("SNP2", "POS2", "CHR2")])
    snp2pos$pos_index2 <- seq.int(nrow(snp2pos))
    snp2merge <- snp2pos[, colnames(snp2pos) %in% c("SNP2", "pos_index2")]
    d_order2 <- merge(d_order1, snp2merge, by="SNP2")
    plotdat <- d_order2
  } else {
    #Duplicate columns and bind with original dataset
    snp1pos <- unique(d_order[, colnames(d_order) %in% c("SNP1", "POS1", "CHR1")])
    snp2pos <- unique(d_order[, colnames(d_order) %in% c("SNP2", "POS2", "CHR2")])
    names(snp1pos) <- c("SNP", "CHR", "POS")
    names(snp2pos) <- c("SNP", "CHR", "POS")
    snppos <- rbind(snp1pos, snp2pos)
    snporder <- unique(snppos[order(snppos$CHR, snppos$POS), ])
    snporder$pos_index <- seq.int(nrow(snporder))
    snpmerge <- snporder[, c(1,4)]
    names(snpmerge) <- c("SNP1", "pos_index1")
    dat1 <- merge(d, snpmerge, by="SNP1", all=TRUE)
    names(snpmerge) <- c("SNP2", "pos_index2")
    dat2 <- merge(dat1, snpmerge, by="SNP2", all=TRUE)
    dat3 <- dat2
    names(dat3)[names(dat3) == 'SNP1'] <- 'placeholder'
    names(dat3)[names(dat3) == 'SNP2'] <- 'SNP1'
    names(dat3)[names(dat3) == 'placeholder'] <- 'SNP2'
    plotdat <- rbind(dat2, dat3)
  }

  #Set tooltip
  xsub <- plotdat[colnames(plotdat) %in% c("SNP1", "CHR1", "POS1","pos_index1")]
  ysub <- plotdat[colnames(plotdat) %in% c("SNP2", "CHR2", "POS2","pos_index2")]

  #Set up dataframe with chromosome position info x axis
  xmaxRows <- by(xsub, xsub$CHR1, function(x) x[which.max(x$pos_index1),])
  xminRows <- by(xsub, xsub$CHR1, function(x) x[which.min(x$pos_index1),])
  xmilimits <- do.call(rbind, xminRows)
  xmalimits <- do.call(rbind, xmaxRows)
  xlims <- merge(xmilimits, xmalimits, by="CHR1")
  names(xlims) <- c("CHR", "snpx", "posx", "posmin", "snpy", "posy", "posmax")
  xlims$av <- (xlims$posmin + xlims$posmax)/2
  xlims <- xlims[order(xlims$CHR),]

  #Set up dataframe with chromosome position info x axis
  ymaxRows <- by(ysub, ysub$CHR2, function(x) x[which.max(x$pos_index2),])
  yminRows <- by(ysub, ysub$CHR2, function(x) x[which.min(x$pos_index2),])
  ymilimits <- do.call(rbind, yminRows)
  ymalimits <- do.call(rbind, ymaxRows)
  ylims <- merge(ymilimits, ymalimits, by="CHR2")
  names(ylims) <- c("CHR", "snpx", "posx", "posmin", "snpy", "posy", "posmax")
  ylims$av <- (ylims$posmin + ylims$posmax)/2
  ylims <- ylims[order(ylims$CHR),]


  if(!is.null(me)){
    me$Frame <- factor(me$Frame)
    #Bind with main effect
    colnames(me)[1] <- "SNP1"
    xme <- unique(merge(me, plotdat[, colnames(plotdat) %in% c("SNP1", "pos_index1")], by="SNP1"))
    colnames(me)[1] <- "SNP2"
    yme <- unique(merge(me, plotdat[, colnames(plotdat) %in% c("SNP2", "pos_index2")], by="SNP2"))
  }

  #Plot
  p <- ggplot(plotdat, aes(x=pos_index1, y=pos_index2, fill=pvalue, frame=Frame)) + geom_tile_interactive()
  p <- p + scale_x_continuous(breaks=xlims$posmax+0.5, labels=xlims$CHR, expand=c(0,0))
  p <- p + scale_y_continuous(breaks=ylims$posmax+0.5, labels=ylims$CHR, expand=c(0,0))
  if(!is.null(me)){
    if("Shape" %in% names(me)){
      xme$Shape <- factor(xme$Shape)
      yme$Shape <- factor(yme$Shape)
      p <- p + geom_point(data=xme, aes(x=pos_index1, y= -1, color=pvalue, shape=Shape, frame=Frame)) + expand_limits(y=-2) + geom_point(data=yme, aes(x=-1, y=pos_index2, color=pvalue, shape=Shape, frame=Frame)) + expand_limits(x=-2)
      p <- p + labs(shape="")
    } else {
      p <- p + geom_point(data=xme, aes(x=pos_index1, y= -1, color=pvalue, frame=Frame)) + expand_limits(y=-2) + geom_point(data=yme, aes(x=-1, y=pos_index2, color=pvalue, frame=Frame)) + expand_limits(x=-2)
    }
    p <- p + scale_color_gradientn(colours=color_vector, values=value_vector, guide = FALSE)
    p <- p + guides(fill = guide_legend(override.aes = list(shape = NA)))
  }

  #Add colorbar or legend
  if(highlight_p!="off"){
    p <- p + scale_fill_gradientn(colours=color_vector, values=value_vector, name=legend,limits=c(0,1),breaks=breaks_vector, labels=label_vector, guide="legend")
  } else {
    p <- p + scale_fill_gradientn(colours=color_vector, values=value_vector, name=legend,limits=c(0,1),breaks=breaks_vector, labels=label_vector, guide="colorbar")
  }

  #Format
  p <- p + theme(panel.grid.minor = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
  p <- p + ggtitle(title)

  #Animate and save
  print(paste("Saving plot to ", file, ".", ext, sep=""))
  ap <- gganimate(p)
  gganimate_save(ap, filename=paste(file, ".", ext, sep=""))
  return(ap)

}
