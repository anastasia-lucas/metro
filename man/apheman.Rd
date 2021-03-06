% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apheman.R
\name{apheman}
\alias{apheman}
\title{apheman}
\usage{
apheman(d, phegroup, line, log10 = TRUE, yaxis, ymax, opacity = 1,
  annotate_snp, annotate_p, highlight_snp, highlight_p, highlighter = "red",
  title = NULL, chrcolor1 = "#AAAAAA", chrcolor2 = "#4D4D4D", groupcolors,
  background = "variegated", chrblocks = TRUE, file = "apheman",
  ext = "gif", hgt = 800, wi = 1300)
}
\arguments{
\item{d}{data frame, must contain PHE, SNP, CHR, POS, pvalue, Frame columns, optional Shape}

\item{phegroup}{optional grouping file for phenotypes, must contain PHE and Group columns}

\item{line}{optional pvalue threshold to draw red line at}

\item{log10}{plot -log10() of pvalue column, boolean}

\item{yaxis}{label for y-axis, automatically set if log10=TRUE}

\item{ymax}{set the upper limit for the y-axis if not automatically scaled}

\item{opacity}{opacity of points, from 0 to 1, useful for dense plots}

\item{annotate_snp}{vector of SNPs to annotate}

\item{annotate_p}{pvalue threshold to annotate}

\item{highlight_snp}{vector of SNPs to highlight}

\item{highlight_p}{pvalue threshold to highlight}

\item{highlighter}{color to highlight}

\item{title}{optional string for plot title}

\item{chrcolor1}{first alternating color for chromosome}

\item{chrcolor2}{second alternating color for chromosome}

\item{groupcolors}{named vector of colors for data in 'Color' column}

\item{background}{variegated or white}

\item{chrblocks}{boolean, turns on x-axis chromosome marker blocks}

\item{file}{file name of saved image}

\item{ext}{file type to save, "gif" or "mp4"}

\item{hgt}{height of plot in pixels}

\item{wi}{width of plot in pixels}
}
\value{
.gif or .mp4 file
}
\description{
Create Manhattan plots for PheWAS
}
\examples{
#We can duplicate the PHE column so that each
#phenotype will be a frame
data(phewas)
phewas$Frame <- phewas$PHE
apheman(phewas, line=0.001, title="PheWAS Example:")
}
\seealso{
\code{\link{pheman}}, \code{\link{ipheman}}, \code{\link{agman}}, \code{\link{aeman}}

Other PheWAS functions: \code{\link{ipheman}},
  \code{\link{pheman}}

Other animated plotting functions: \code{\link{aeman}},
  \code{\link{agman}}
}
