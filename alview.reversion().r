#last edited at 20191226

#a function with the reverse funtion of ape::alview() and seqinr::dotPlot()
#i.e. interpret dots in 2nd and latter sequences according to the 1st sequence

alview.reversion=function(strings){
  #reasons for why not using ape::DNAbin object to manipulate data, but use additional
  #Biostrings::BStringSet object to manipulate and then convert BStringSet into DNAbin:
  #DNAbin can't store "." in a sequence, however there is doomed to contain "."
  #
  #>>>r package dependency, function dependency, internal functions begin>>>
  library(Biostrings)
  library(magrittr)
  library(stringr)
  library(ape)
  #<<<r package dependency, function dependency, internal functions end<<<
  #>>>control the input of strings begin>>>
  if(!is.vector(strings)|is.list(strings)){stop("strings must be a character vector.")}
  if(length(strings)==0){stop("length(strings) can't be 0.")}
  strings=as.character(strings)
  if(nchar(strings) %>% unique() %>% length()!=1){stop("all nchar(strings) must be the same.")}
  #
  ans=grep(pattern="^ +[[:digit:]]+$",x=strings)
  if(length(ans)>0)if(!all(ans==1:length(ans))){
    stop("strings begin with spaces must represent residue indices. and these strings must occur at first.")
  }
  strings=strings[-ans]
  rm(ans)
  strings=str_match(string=strings,pattern="^([^ ]+) +([^ ]+)$")[,-1]
  if(anyNA(strings[,1])){stop("strings can't match on the valid pattern of result of ape::alview().")}
  if(nchar(strings[,2]) %>% unique() %>% length()!=1){stop("sequence lengths in strings must be same.")}
  strings=setNames(BStringSet(toupper(strings[,2])),strings[,1])#stop when strings contain invalid characters
  #
  l=c(length(strings),nchar(strings)[1])
  if(countPattern(pattern=".",strings[[1]])!=0){stop("characters in the 1st sequence can't be \".\"")}
  #<<<control the input of strings end<<<
  #>>>main manipulation begin>>>
  if(l[1]>1)for(i in 1:l[2]){
    #l[1]>1 controls not entering the loop if strings only contain one sequence
    tulip=IRanges(start=i,width=1) %>% IRangesList() %>% rep(times=l[1])
    dahlia=extractAt(x=strings,at=tulip) %>% as.list()
    for(j in 1:l[1]){dahlia[[j]]=dahlia[[j]][[1]] %>% as.character();}
    for(j in 2:l[1]){
      if(dahlia[[j]]==dahlia[[1]]){
        stop("characters in other sequences (use \".\" if same) can't be same with the 1st sequence.")
      }
      if(dahlia[[j]]=="."){dahlia[[j]]=dahlia[[1]];}
    }
    for(j in 1:l[1]){dahlia[[j]]=BStringSet(dahlia[[j]]);}
    dahlia=BStringSetList(dahlia)
    strings=replaceAt(strings,at=tulip,value=dahlia)
    rm(tulip,dahlia,j)
  }
  #<<<main manipulation end<<<
  #>>>change BStringSet object into DNAbin or AAbin object begin>>>
  tulip=uniqueLetters(x=strings)
  if(all(tulip %in% DNA_ALPHABET)){
    strings=DNAStringSet(strings) %>% ape::as.DNAbin()
  }else{
    if(all(tulip %in% AA_ALPHABET)){
      strings=AAStringSet(strings) %>% ape::as.AAbin()
    }else{stop("strings contain invalid characters, both for DNA and AA.");}
  }
  rm(tulip)
  #<<<change BStringSet object into DNAbin or AAbin object end<<<
  return(strings)#i=0;rm(strings,l,i)
}

if(FALSE){
  #here are the test codes
  strings=readLines("~/../OneDrive/r-zm/testalview.txt")
  alview.reversion(strings)
}
