#last edited at 20200101

#a function with the reverse funtion of ape::alview() and seqinr::dotPlot()
#i.e. interpret dots in 2nd and latter sequences according to the 1st sequence

library(Biostrings)
library(magrittr)
library(stringr)
library(ape)
library(seqinr)
source(file="https://github.com/ywd5/r-zm/raw/master/t.XStringSet().r")

alview.reversion=function(strings){
  #reasons for why not using ape::DNAbin object to manipulate data, but use additional
  #Biostrings::BStringSet object to manipulate and then convert BStringSet into DNAbin:
  #DNAbin can't store "." in a sequence, however there is doomed to contain "."
  #
  #>>>internal functions begin>>>
  #<<<internal functions end<<<
  metainfo=list(type="DNA",nrow=0,ncol=0,fString=DNAString,fStringSet=DNAStringSet)
  #>>>control the input of strings and get metainfo begin>>>
  if(!is.vector(strings)|is.list(strings)){stop("strings must be a character vector.")}
  if(length(strings)==0){stop("length(strings) can't be 0.")}
  tempo=BStringSet(strings[1]);if(length(strings)>1)for(i in 2:length(strings)){
    tempo=c(tempo,BStringSet(strings[i]))
  };strings=tempo;rm(tempo)
  if(nchar(strings) %>% unique() %>% length()!=1){stop("all nchar(strings) must be the same.")}
  #
  #ans=vcountPattern(pattern="^ +[[:digit:]]+$",subject=strings,fixed=FALSE) %>% {(1:length(.))[.==1];}
  ans=rep(0,times=length(strings));for(i in 1:length(strings)){
    ans[i]=str_count(string=as.character(strings[[i]]),pattern="^ +[[:digit:]]+$")
  };ans=(1:length(ans))[ans==1];
  if(length(ans)>0)if(!all(ans==1:length(ans))){
    stop("strings begin with spaces must represent residue indices. and these strings must occur at first.")
  }
  if(length(ans)>0){strings=strings[-ans];}
  rm(ans)
  metainfo$nrow=length(strings)
  #
  tempo=rep("",times=metainfo$nrow)
  for(i in 1:metainfo$nrow){
    ans=as.character(strings[[i]]) %>% str_match(pattern="^([^ ]+ +)([^ ]+)$") %>% .[,-1]
    if(is.na(ans[1])){stop("strings can't match on the valid pattern of result of ape::alview().")}
    tempo[i]=ans[1];strings[[i]]=BString(toupper(ans[2]))
    rm(ans)
  }
  if(nchar(tempo) %>% unique() %>% length()!=1){stop("sequence lengths in strings must be same.")}
  names(strings)=str_trim(string=tempo,side="right")
  rm(tempo)
  metainfo$ncol=nchar(strings[[1]])
  #
  if(countPattern(pattern=".",strings[[1]],fixed=TRUE)!=0){
    stop("characters in the 1st sequence can't be \".\"")
  }
  ans=uniqueLetters(x=strings)
  if(!all(ans %in% DNA_ALPHABET)){
    if(all(ans %in% AA_ALPHABET)){
      metainfo$type="AA";metainfo$fString=AAString;metainfo$fStringSet=AAStringSet;
    }else{
      stop("strings contain invalid characters, both for DNA and AA.")
    }
  }
  strings=metainfo$fStringSet(strings)
  rm(ans)
  #<<<control the input of strings and get metainfo end<<<
  #>>>main manipulation begin>>>
  if(metainfo$nrow>1){
    tulip=t.XStringSet(xss=strings)
    for(i in 1:length(tulip)){
      dahlia=as.character(tulip[[i]])
      if(str_count(string=dahlia,pattern=substr(dahlia,1,1))>1){
        stop("characters in other sequences (use \".\" if same) can't be same with the 1st sequence.")
      }
      tulip[[i]]=chartr(old=".",new=substr(dahlia,1,1),x=dahlia) %>% metainfo$fString()
      rm(dahlia)
    }
    strings=t.XStringSet(xss=tulip) %>% setNames(object=.,nm=names(strings))
    rm(tulip)
  }
  #if(metainfo$nrow>1)for(i in 1:metainfo$ncol){
  #  #metainfo$nrow>1 controls not entering the loop if strings only contain one sequence
  #  tulip=IRanges(start=i,width=1) %>% IRangesList() %>% rep(times=metainfo$nrow)
  #  dahlia=extractAt(x=strings,at=tulip) %>% as.list()
  #  for(j in 1:metainfo$nrow){dahlia[[j]]=dahlia[[j]][[1]] %>% as.character();}
  #  for(j in 2:metainfo$nrow){
  #    if(dahlia[[j]]==dahlia[[1]] & dahlia[[1]]!="-"){
  #      stop("characters in other sequences (use \".\" if same) can't be same with the 1st sequence.")
  #    }
  #    if(dahlia[[j]]=="."){dahlia[[j]]=dahlia[[1]];}
  #  }
  #  for(j in 1:metainfo$nrow){dahlia[[j]]=BStringSet(dahlia[[j]]);}
  #  dahlia=BStringSetList(dahlia)
  #  strings=replaceAt(strings,at=tulip,value=dahlia)
  #  rm(tulip,dahlia,j)
  #}
  #<<<main manipulation end<<<
  if(metainfo$type=="DNA"){strings=ape::as.DNAbin(strings);}else{strings=ape::as.AAbin(strings);}
  strings=as.matrix(strings)
  i=0;rm(metainfo,i)
  return(strings)
}

if(FALSE){
  #here are the test codes
  strings=c("   1234","x  acgt","y  ...g","xy ..c.")
  alview.reversion(strings)
}
