#last edited at 20200102

#a function to convert ape::DNAbin object to Biostrings::DNAStringSet object

library(Biostrings)
library(ape)
library(magrittr)

as.DNAStringSet.DNAbin <- function(dnabin){
  #>>>control the input of dnabin begin>>>
  if(!is(object=dnabin,class2="DNAbin")){stop("dnabin must be an DNAbin object defined by ape.");}
  if(is.matrix(dnabin)){dnabin=as.list(dnabin);}
  if(!is.list(dnabin)){stop("dnabin must be in list or matrix form, but not vector or other form.")}
  if(length(dnabin)==0){stop("dnabin must contain at least 1 sequence.")}
  #<<<control the input of dnabin end<<<
  #>>>main manipulation begin>>>
  l=length(dnabin)
  result=DNAStringSet("") %>% rep(times=l) %>% setNames(nm=labels(dnabin))
  for(i in 1:l)if(anyNA(as.character(dnabin[i])[[1]])){stop("NO. ",i,"-th sequence in dnabin contain NA.")}
  for(i in 1:l){result[i]=as.character(dnabin[i]) %>% .[[1]] %>% paste0(collapse="")}
  #<<<main manipulation end<<<
  rm(dnabin,l,i)
  return(result)
}

if(FALSE){
  #here are the test codes
  dnabin=list(c("a","c","g","t","M","r"),"") %>% as.DNAbin()
  dnabin=list(c("a","c","g","t","M","r"),character()) %>% as.DNAbin()
  dnabin=list(c("a","c","g","t","M","r"),character()) %>% as.DNAbin() %>% setNames(nm=c("x","y"))
  as.DNAStringSet(dnabin)
}
