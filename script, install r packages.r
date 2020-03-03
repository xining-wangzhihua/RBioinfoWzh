#last edited at 20200303


message("before installing r packages, make sure you have installed Rtools and",
        "jre (java runtime environment) manually, and set path to Rtools in r.")
#
#try this code to set the path to Rtools:
#Sys.setenv(BINPREF="D:/softwares/Rtools/mingw_$(WIN)/bin/")
#
update.packages()


l.common=c("devtools","ggplot2","gtools","Hmisc","magrittr","questionr",
           "R.methodsS3","R.oo","R6","Rcpp","R.utils","readxl","rlang","rJava","shiny",
           "stringr","tibble","xml2")
#Hmisc may be installed only for `%nin%`
l.rare=c("curl","data.table","digest","dplyr","fastmap","fastmatch","gdata",
         "measurements","OSMscale","plyr","RCurl","readr","rlist","rncl",
         "RNeXML","stringi","tidyr","utf8","vctrs","XiMpLe","xlsx","XML")
l.bioconductor=c("Biobase","BiocGenerics","IRanges","GenomicRanges","S4Vectors",
                 "XVector","Biostrings","sangerseqR","msa")
l.phylo.and.bioinfo=c("ape","haplo.stats","phylobase","seqinr")#doesn't include packages in bioconductor
l.geography=c("mapdata","mapproj","maps","maptools","OpenStreetMap","rgeos","sf")


l.common           =setdiff(x=l.common,           y=rownames(installed.packages()))
l.phylo.and.bioinfo=setdiff(x=l.phylo.and.bioinfo,y=rownames(installed.packages()))
l.geography        =setdiff(x=l.geography,        y=rownames(installed.packages()))
if(length(c( l.common, l.phylo.and.bioinfo, l.geography ))>0){
  install.packages(pkgs=c( l.common, l.phylo.and.bioinfo, l.geography ))
}
#
if( ! "BiocManager" %in% rownames(installed.packages()) ){
  install.packages(pkgs="BiocManager")
}
l.bioconductor=setdiff(x=l.bioconductor,y=rownames(installed.packages()))
if(length(l.bioconductor)>0){BiocManager::install(pkgs=l.bioconductor);}
#
if(FALSE){
  l.rare=setdiff(x=l.rare,y=rownames(installed.packages()))
  if(length(l.rare)>0){install.packages(pkgs=l.rare);}
}


rm(l.common, l.rare, l.bioconductor, l.phylo.and.bioinfo, l.geography)
