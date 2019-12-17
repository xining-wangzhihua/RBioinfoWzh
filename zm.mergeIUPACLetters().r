#last edited at 20191217

#a function similar to Biostrings::mergeIUPACLetters, but accept 2 strings or DNAString object as 2 sequences.

zm.mergeIUPACLetters=function(seq1,seq2){
  #>>>r package dependency, function dependency, internal functions begin>>>
  library(Biostrings)
  library(tibble)
  library(magrittr)
  #<<<r package dependency, function dependency, internal functions end<<<
  #>>>control the input of seq1 and seq2 begin>>>
  if(is(object=seq1,class2="DNAString")){seq1=as.character(seq1);}
  if(is(object=seq2,class2="DNAString")){seq2=as.character(seq2);}
  if(!is.vector(seq1)|is.list(seq1)|!is.vector(seq2)|is.list(seq2)){
    stop("seq1 and seq2 must all be character string or DNAString objects.")
  }
  seq1=as.character(seq1)
  seq2=as.character(seq2)
  if(length(seq1)!=1|length(seq2)!=1){stop("seq1 and seq2 must all contain 1 sequence.");}
  if(nchar(seq1)!=nchar(seq2)|seq1==""){stop("residue numbers in seq1 and seq2 must equal, and can't be 0.");}
  #<<<control the input of seq1 and seq2 end<<<
  #>>>main manipulation begin>>>
  seq1=strsplit(seq1,split="")[[1]] %>% toupper()
  seq2=strsplit(seq2,split="")[[1]] %>% toupper()
  if(unique(seq1,seq2) %>% setdiff(x=.,y=DNA_ALPHABET) %>% length()>0){
    stop("sequences in seq1 or seq2 contain invalid characters for DNA.")
  }
  for(i in 1:length(seq1))if(seq1[i]!=seq2[i]){
    msu=seq1[i];spbu=seq2[i];
    if((msu=="-"|msu=="+"|msu==".") & (spbu=="-"|spbu=="+"|spbu==".")){
      stop("when residues in seq1 and seq2 are [-+.], they must be the same.")
    }
    if((msu=="-"|msu=="+"|msu==".") & !(spbu=="-"|spbu=="+"|spbu==".")){seq1[i]=seq2[i];}
    if(!(msu=="-"|msu=="+"|msu==".") & (spbu=="-"|spbu=="+"|spbu==".")){}#do nothing
    if(!(msu=="-"|msu=="+"|msu==".") & !(spbu=="-"|spbu=="+"|spbu==".")){
      seq1[i]=paste0(seq1[i],seq2[i]) %>% mergeIUPACLetters()
    }
    rm(msu,spbu)
  }
  seq1=paste0(seq1,collapse="")
  #<<<main manipulation end<<<
  return(seq1)#rm(seq1,seq2,i)
}

if(FALSE){
  #here are the test codes
  ans="";tempo="";
  ans="acg";tempo="acgt";
  ans="acgt";tempo="tgca";
  ans=Biostrings::DNAString("a-c+g.t");tempo="-ac-a.t";
  ans=Biostrings::DNAString("a-c+g.t");tempo="-ac+a.t";
  zm.mergeIUPACLetters(seq1=ans,seq2=tempo)
}
