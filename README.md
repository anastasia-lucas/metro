# metro
An R package for creating Manhattan plots for 'ome-wide association studies

## Overview
metro is an R package for creating static, interactive, and animated Manhattan plots. The package includes functions to visualize genome-wide, phenome-wide, and environment-wide association analysis (GWAS, PheWAS, EWAS, respectively) results, though they may adaptable for other types of data such as Beta, SNP intensity value, or even other types of analyses. This packages allows the user to add meta information which helps the user and collaboraters to better interpret the results.

## Installation
As of now, there is only a development version of the package which can be installed using devtools.

```devtools::install_github('anastasia-lucas/metro')```

This package requires ggplot2, ggiraph for interactive plots, and gganimate (which requires ImageMagick) for animated plots. I recommend Cactuslab's ImageMagick installer, which can be found at http://cactuslab.com/imagemagick/ if you are having issues installing on MacOS. ggrepel is suggested for improved text annotation, but not required.

## Usage

### Create an animated Manhattan plot using GWAS data

```
library(metro)
#Generate some data
#To use the animated plot function, we need to add an animation 'Frame' column to our data
#In this case we will imagine that we've run a GWAS using additive, dominant, and recessive models
#and want to highlight a SNP of interest to see how the p-value changes
add <- data.frame(SNP=paste("rs", seq(1:5000), sep=""),
                  CHR=rep(c(1:22, "X", "Y"), length.out=5000, each=200),
                  POS=rep(seq(1, 10000, by = 200), length.out=5000),
                  pvalue=runif(n=5000),
                  Frame="Additive")
dom <- data.frame(SNP=paste("rs", seq(1:5000), sep=""),
                  CHR=rep(c(1:22, "X", "Y"), length.out=5000, each=200),
                  POS=rep(seq(1, 10000, by = 200), length.out=5000),
                  pvalue=runif(n=5000),
                  Frame="Dominant")
rec <- data.frame(SNP=paste("rs", seq(1:5000), sep=""),
                  CHR=rep(c(1:22, "X", "Y"), length.out=5000, each=200),
                  POS=rep(seq(1, 10000, by = 200), length.out=5000),
                  pvalue=runif(n=5000),
                  Frame="Recessive")
dat <- rbind(add, dom, rec)

agman(d=dat, line=0.0005, highlight_snp="rs4204", chrcolor1="#D4CAA0", chrcolor2="#B3BC92", highlighter="black",opacity=0.7, wi=750, hgt=500)
```
![Alt text](https://media.giphy.com/media/6CBiR3JABocVfJx8pO/giphy.gif)

### Create an interactive Manhattan plot using PheWAS data

```
library(metro)
#Generate some data
#In this case we'd like to also see the p-value when we hover over a point, so we'll add an 'Info' column to the data
#We'd also like to search dbSNP when we click on a point
dat <- data.frame(PHE=rep(paste("Pheno", seq(1:5), sep="") , length.out=5000, each=1),
                  SNP=paste("rs", seq(1:5000), sep=""),
                  CHR=rep(c(1:22, "X", "Y"), length.out=5000, each=200),
                  POS=rep(seq(1, 10000, by = 200), length.out=5000),
                  pvalue=runif(n=5000))
dat$Info <- paste0("p-value:", signif(dat$pvalue, digits=3))

ipheman(d=dat, moreinfo = TRUE, db="dbSNP", line=0.001)
```
![Alt text](https://media.giphy.com/media/ujvtTHsrnXeo5wuZDn/giphy.gif)

### Create a basic Manhattan plot using EWAS data

```
library(metro)
#Generate some data
dat <- data.frame(Variable=paste("Var", seq(1:5000), sep=""), 
                  pvalue=runif(n=5000), 
                  Group=rep(paste("G", seq(1:6), sep=""), length.out=5000, each=1),
                  Shape=rep(paste("S", seq(1:5), sep="") , length.out=5000, each=1))
                  
eman(d=dat, title="EWAS", line=0.001, annotate_p=0.001, color1="#A23B72", color2="#2A84AA", highlight_p=0.001, highlighter="green")
```

![Imgur](https://i.imgur.com/eSARizX.jpg)
