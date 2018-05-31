# metro
An R package for creating Manhattan plots for 'ome-wide association studies

## Overview
metro is an R package for creating static, interactive, and animated Manhattan plots. The package includes functions to visualize genome-wide, phenome-wide, and environment-wide association analysis (GWAS, PheWAS, EWAS, respectively) results, though they may adaptable for other types of data such as Beta, SNP intensity value, or even other types of analyses. THis packages allows the user to add meta information which helps the user and collaboraters to better interpret the results.

## Installation
As of now, there is only a development version of the package which can be installed using devtools.

```devtools::install_github('anastasia-lucas/metro')```

This package requires ggplot2, ggiraph for interactive plots, and gganimate (which requires ImageMagick) for animated plots. I recommend Cactuslab's ImageMagick installer, which can be found at http://cactuslab.com/imagemagick/ if you are having issues installing on MacOS. ggrepel is suggested for improved text annotation, but not required.

## Usage
