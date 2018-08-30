#' iqqunif
#'
#' Create interactive qqplots with an assumed uniform distribution
#' @param d dataframe with at least two columns, p-value and Name; Color, Shape, Info optional
#' @param CI two-sided confidence interval, default 0.95
#' @param splitby if data contains Color and/or Shape, indicate variable(s) by which the data should be subsetted for calculating CIs
#' @param opacity point opacity, default 1
#' @param title plot title
#' @param db choose database to connect to ("dbSNP", "GWASCatalog", or enter your own search address)
#' @param moreinfo includes more information on hover, refers to Info column
#' @param groupcolors named vector of colors corresponding to data in Group column
#' @param highlight_name vector of names to highlight, dataframe must include a Name column
#' @param highlight_p p-value threshold to highlight
#' @param highlighter highlighter color
#' @param annotate_name vector of names to annotate, dataframe must include a Name column
#' @param annotate_p p-value threshold to annotate, dataframe must include a Name column
#' @param line draw a red line at pvalue threshold (observed)
#' @param background can change to "white"
#' @param bigrender can set to TRUE for big plots (~50000 rows) that produce huge input lookup error
#' @param file filename
#' @param wi width of plot
#' @param hgt height of plot
#' @param res resolution of plot
#' @return html file
#' @import ggplot2
#' @export
#' @examples
#' data(gwas)
#' qqdat <- data.frame(pvalue=gwas$pvalue, Color=gwas$Frame, Name=gwas$Frame)
#' iqqunif(d=qqdat, splitby="Color", db="https://www.google.com/search?q=")

iqqunif <- function(d, CI=0.95, opacity=1, groupcolors, splitby=NULL, moreinfo=TRUE, db, highlight_p, highlight_name, annotate_p, annotate_name, highlighter="red", line, background, title, bigrender=FALSE, file="iqqunif", wi=7, hgt=7, res=300){
  if (!requireNamespace(c("ggiraph"), quietly = TRUE)==TRUE) {
    stop("Please install ggiraph to create interactive visualizations.", call. = FALSE)
  }

  if("Color" %in% colnames(d)){
    if(!missing(groupcolors)){
      colrs <- groupcolors
    } else {
      ngroupcolors <- nlevels(factor(d$Color))
      if(ngroupcolors > 15){
        if (!requireNamespace(c("RColorBrewer"), quietly = TRUE)==TRUE) {
          stop("Please install RColorBrewer to add color attribute for more than 15 colors.", call. = FALSE)
        } else {
          getPalette = grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))
          colrs<- getPalette(ngroupcolors)
        }
      } else {
        pal <- pal <- c("#009292", "#920000", "#490092", "#db6d00", "#24ff24",
                        "#ffff6d", "#000000", "#006ddb", "#004949","#924900",
                        "#ff6db6", "#6db6ff","#b66dff", "#ffb6db","#b6dbff")
        colrs <- pal[1:ngroupcolors]
      }
    }
  }

  if(!is.null(splitby)){
    dlist <- split(d, d[, splitby])
    df <- lapply(dlist, function(x) cbind(x[order(x$pvalue),],
                                          obs=-log10(sort(x$pvalue)),
                                          ex=-log10(stats::ppoints(length(!is.na(x$pvalue)))),
                                          cl=-log10(stats::qbeta(p = (1-CI)/2, shape1 = 1:length(!is.na(x$pvalue)), shape2 = length(!is.na(x$pvalue)):1)),
                                          cu=-log10(stats::qbeta(p = (1+CI)/2, shape1 = 1:length(!is.na(x$pvalue)), shape2 = length(!is.na(x$pvalue)):1))))
    dat <- do.call("rbind", df)
  } else {
    dat <- cbind(d[order(d$pvalue), , drop=FALSE],
                 obs=-log10(sort(d$pvalue)),
                 ex=-log10(stats::ppoints(length(!is.na(d$pvalue)))),
                 cl=-log10(stats::qbeta(p = (1-CI)/2, shape1 = 1:length(!is.na(d$pvalue)), shape2 = length(!is.na(d$pvalue)):1)),
                 cu=-log10(stats::qbeta(p = (1+CI)/2, shape1 = 1:length(!is.na(d$pvalue)), shape2 = length(!is.na(d$pvalue)):1)))
  }

  #Set up tooltip
  ###See what info would be useful here i.e. SNP or something else
  dat$tooltip <- if (moreinfo==TRUE) c(paste0(dat$Name, "\n ", dat$Info)) else dat$Name

  #Set up onclick
  if(!missing(db)){
    if(db=="GWASCatalog"){
      dat$onclick <- sprintf("window.open(\"%s%s\")","https://www.ebi.ac.uk/gwas/search?query=", as.character(dat$Name))
    } else if(db=="dbSNP"){
      dat$onclick <- sprintf("window.open(\"%s%s\")","https://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?searchType=adhoc_search&type=rs&rs=", as.character(dat$Name))
    } else {
      dat$onclick <- sprintf("window.open(\"%s%s\")", db, as.character(dat$Name))
    }
  } else {
    dat$onclick <- NA
  }

  if("Shape" %in% splitby & "Color" %in% splitby){
    linaes1 <- "geom_line(aes(ex, cu, linetype=Shape, color=Color))"
    linaes2 <- "geom_line(aes(ex, cl, linetype=Shape, color=Color))"
  } else if("Shape" %in% splitby) {
    linaes1 <- "geom_line(aes(ex, cu, linetype=Shape))"
    linaes2 <- "geom_line(aes(ex, cl, linetype=Shape))"
  } else if("Color" %in% splitby){
    linaes1 <- "geom_line(aes(ex, cu, color=Color))"
    linaes2 <- "geom_line(aes(ex, cl, color=Color))"
  } else {
    linaes1 <- "geom_line(aes(ex, cu), linetype=2)"
    linaes2 <- "geom_line(aes(ex, cl), linetype=2)"
  }
  #Plot
  if("Shape" %in% colnames(d)){
    if("Color" %in% colnames(d)){
      p <- ggplot(dat, aes(ex, obs)) + ggiraph::geom_point_interactive(aes(tooltip=tooltip, onclick=onclick, shape=Shape, colour=Color),alpha=opacity)
    } else {
      p <- ggplot(dat, aes(ex, obs)) + ggiraph::geom_point_interactive(aes(tooltip=tooltip, onclick=onclick, shape=Shape),alpha=opacity)
    }
    p <- p + theme(legend.position="bottom", legend.title = element_blank())
  } else {
    if("Color" %in% colnames(d)){
      p <- ggplot(dat, aes(ex, obs)) + ggiraph::geom_point_interactive(aes(tooltip=tooltip, onclick=onclick, colour=Color), alpha=opacity)
      p <- p + theme(legend.position="bottom", legend.title = element_blank())
    } else{
      p <- ggplot(dat, aes(ex, obs)) + ggiraph::geom_point_interactive(aes(tooltip=tooltip, onclick=onclick), alpha=opacity)
    }
  }
  p <- p + eval(parse(text=linaes1)) + eval(parse(text=linaes2))
  p <- p + labs(x=expression(paste("Expected -log"[10], plain(P))), y=expression(paste("Observed -log"[10], plain(P))))
  p <- p + geom_abline(intercept = 0, slope = 1, alpha = 0.5)
  p <- p + theme(panel.grid.minor = element_blank())
  if("Color" %in% colnames(d)){p <- p + scale_color_manual(name = "Color", values = colrs)}
  if(!missing(line)){p <- p + geom_hline(yintercept = redline, colour="red")}
  if(!missing(title)) {p <- p + ggtitle(title)}
  if(!missing(background)) {p <- p + theme(panel.background = element_rect(fill=background))}
  #Add extra annotations
  if(!missing(highlight_name)){
    if("Shape" %in% names(d)){
      p <- p + ggiraph::geom_point_interactive(data=dat[dat$Name %in% highlight_name, ], aes(x=ex, y=obs, tooltip=tooltip, onclick=onclick, shape=Shape), colour=highlighter)
      p <- p + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p <- p + ggiraph::geom_point_interactive(data=dat[dat$Name %in% highlight_name, ], aes(x=ex, y=obs, tooltip=tooltip, onclick=onclick), colour=highlighter)
    }
  }
  if(!missing(highlight_p)){
    if("Shape" %in% names(d)){
      p <- p + ggiraph::geom_point_interactive(data=dat[dat$pvalue < highlight_p, ], aes(x=ex, y=obs, tooltip=tooltip, onclick=onclick, shape=Shape), colour=highlighter)
      p <- p + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p <- p + ggiraph::geom_point_interactive(data=dat[dat$pvalue < highlight_p, ], aes(x=ex, y=obs, tooltip=tooltip, onclick=onclick), colour=highlighter)
    }
  }
  if(!missing(annotate_p)){
    if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE) {
      print("Consider installing 'ggrepel' for improved text annotation")
      p <- p + geom_text(data=dat[dat$pvalue < annotate_p,], aes(ex,obs,label=Name), color="black")
    } else {
      p <- p + ggrepel::geom_text_repel(data=dat[dat$pvalue < annotate_p,], aes(ex,obs,label=Name), color="black")
    }
  }
  if(!missing(annotate_name)){
    if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE){
      print("Consider installing 'ggrepel' for improved text annotation")
      p <- p + geom_text(data=dat[dat$Name %in% annotate_name,], aes(ex,obs,label=Name))
    } else {
      p <- p + ggrepel::geom_text_repel(data=dat[dat$Name %in% annotate_name,], aes(ex,obs,label=Name))
    }
  }
  #Save
  print(paste("Saving plot to ", file, ".html", sep=""))
  tooltip_css <- "background-color:black;color:white;padding:6px;border-radius:15px 15px 15px 15px;"
  if(bigrender==TRUE){
    print(paste("WARNING: Attempting to render ", nrow(dat), " rows. Plot may be slow to render or load.", sep=""))
    ip <- ggiraph::ggiraph(code=print(p), tooltip_extra_css = tooltip_css, tooltip_opacity = 0.75, zoom_max = 6, width_svg=wi, height_svg=hgt, xml_reader_options = list(options="HUGE"))
    htmlwidgets::saveWidget(widget=ip, file=paste(file, ".html", sep=""))
    return(p)
  } else {
    ip <- ggiraph::ggiraph(code=print(p), tooltip_extra_css = tooltip_css, tooltip_opacity = 0.75, zoom_max = 6, width_svg=wi, height_svg=hgt)
    htmlwidgets::saveWidget(widget=ip, file=paste(file, ".html", sep=""))
    return(ip)
  }
}
