% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gman.R
\name{gman}
\alias{gman}
\title{gman}
\usage{
gman(d, line, log10 = TRUE, yaxis, ymax, opacity = 1, annotate_snp,
  annotate_p, highlight_snp, highlight_p, highlighter = "red", title = NULL,
  chrcolor1 = "#AAAAAA", chrcolor2 = "#4D4D4D", groupcolors,
  background = "variegated", chrblocks = FALSE, file = "gman", hgt = 7,
  wi = 12, res = 300)
}
\arguments{
\item{d}{data frame, must contain SNP, CHR, POS, pvalue columns, optional Shape and Color}

\item{line}{optional pvalue threshold to draw red line at}

\item{log10}{plot -log10() of pvalue column, boolean}

\item{yaxis}{label for y-axis, automatically set if log10=TRUE}

\item{ymax}{set the upper limit for the y-axis if not automatically scaled}

\item{opacity}{opacity of points, from 0 to 1, useful for dense plots}

\item{annotate_snp}{vector of SNPs to annotate}

\item{annotate_p}{pvalue threshold to annotate}

\item{highlight_snp}{vector of snps to highlight}

\item{highlight_p}{pvalue threshold to highlight}

\item{highlighter}{color to highlight}

\item{title}{optional string for plot title}

\item{chrcolor1}{first alternating color for chromosome}

\item{chrcolor2}{second alternating color for chromosome}

\item{groupcolors}{named vector of colors for data in 'Color' column}

\item{background}{variegated or white}

\item{chrblocks}{boolean, turns on x-axis chromosome marker blocks}

\item{file}{file name of saved image}

\item{hgt}{height of plot in inches}

\item{wi}{width of plot in inches}

\item{res}{resolution of plot in pixels per inch}
}
\value{
png image
}
\description{
Create Manhattan plots for GWAS
}
\examples{
data(gwas)
gman(d=gwas[gwas$Frame=="Additive", 1:4], line=0.0005, title="GWAS Example: Additive")
}
\seealso{
\code{\link{igman}}, \code{\link{agman}}, \code{\link{pheman}}, \code{\link{eman}}

Other GWAS functions: \code{\link{agman}},
  \code{\link{igman}}

Other static plotting functions: \code{\link{eman}},
  \code{\link{pheman}}
}
