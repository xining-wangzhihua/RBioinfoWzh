#last edited at 20191217

#a function with the reverse funtion of ape::alview()
#i.e. interpret dots in the target sequence according to a reference sequence

reverse.alview=function(seq,ref){
  #>>>r package dependency, function dependency, internal functions begin>>>
  library(Biostrings)
  #<<<r package dependency, function dependency, internal functions end<<<
  #>>>control the input of seq and ref begin>>>
  if(is(object=seq,class2="XString")){seq=as.character(seq);}
  if(is(object=ref,class2="XString")){ref=as.character(ref);}
  if(!is.vector(seq)|is.list(seq)|!is.vector(ref)|is.list(ref)){
    stop("seq and ref must be character string or DNAString objects.")
  }
  if(length(seq)!=1|length(ref)!=1){stop("seq and ref must all contain one sequence.");}
  if(nchar(seq)!=nchar(ref)|seq==""){stop("residue numbers in seq and ref must equal, and can't be 0.");}
  #<<<control the input of seq and ref end<<<
  #>>>main manipulation begin>>>
  seq=strsplit(seq,split="")[[1]] %>% toupper()
  ref=strsplit(ref,split="")[[1]] %>% toupper()
  for(i in 1:length(seq))if(seq[i]=="."){seq[i]=ref[i];}
  seq=paste0(seq,collapse="")
  #<<<main manipulation end<<<
  return(seq)#rm(seq,ref,i)
}

if(FALSE){
  #here are the test codes
  reverse.alview(seq="a.b.c",ref="-;/!~")
}
