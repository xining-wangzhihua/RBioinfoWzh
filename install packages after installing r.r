#last edited at 20191211

list_1=c("ape","dplyr","fastmatch","ggplot2","gtools","magrittr","questionr",
         "readr","rlang","shiny","stringi","stringr","tibble","XML","XiMpLe")
list_bioconductor=c("Biostrings","msa","S4Vectors","sangerseqR")
list_2=c("curl","data.table","devtools","digest","fastmap","plyr","R.oo","seqinr",
         "tidyr","utf8","vctrs")

install.packages(pkgs=list_1)
install.packages(pkgs="BiocManager")
BiocManager::install(pkgs=list_bioconductor)
if(FALSE){install.packages(pkgs=list_2);}
