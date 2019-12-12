#last edited at 20191211

update.packages()

list_1=c("ape","dplyr","fastmatch","ggplot2","gtools","magrittr","questionr",
         "readr","rlang","shiny","stringi","stringr","tibble","XML")
list_bioconductor=c("Biostrings","msa","S4Vectors","sangerseqR")
list_2=c("curl","data.table","devtools","digest","fastmap","plyr","R.oo","seqinr",
         "tidyr","utf8","vctrs","XiMpLe")

install.packages(pkgs=list_1)
install.packages(pkgs="BiocManager")
BiocManager::install(pkgs=list_bioconductor)
if(FALSE){install.packages(pkgs=list_2);}

rm(list_1,list_2,list_bioconductor)
