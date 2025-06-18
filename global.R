.libPaths(c("/usr/lib/R/site-library",.libPaths()))
options(java.parameters = "-Xmx4g")
library(rJava)
.jinit() 
library(rTASSEL) #
library(shinycssloaders)
library(VariantTools) #
library(Rsamtools) #
library(rtracklayer) #
library(gmapR) #
library(ShortRead) #
library(Biostrings)
library(stringr)
library(rvest)
library(ggplot2) #
library(qtl) #
library(ASMap) #
library(gridExtra) #for stacking the plots
library(parallel)
library(qtlhot) #
library(LinkageMapView) #
library(vcfR) #for vcf plotting
#library(tidyr) #for vcf plotting
library(dplyr) #for plotting vcf
library(zip)
#library(shinyFiles)
#library(viridis) #For color palette
library(ggrepel) #For labelling the marker
library(reshape2)
