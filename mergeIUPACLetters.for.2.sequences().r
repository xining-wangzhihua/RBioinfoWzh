#last edited at 20191227

#a function similar to Biostrings::mergeIUPACLetters, but accept 2 strings representing 2 sequences.

mergeIUPACLetters.for.2.sequences=function(seq1,seq2){
  #>>>r package dependency, function dependency, internal functions begin>>>
  library(Biostrings)
  library(magrittr)
  #<<<r package dependency, function dependency, internal functions end<<<
  #>>>control the input of seq1 and seq2 begin>>>
  if(!is.vector(seq1)|is.list(seq1)|!is.vector(seq2)|is.list(seq2)){
    stop("seq1 and seq2 must all be character strings.")
  }
  seq1=as.character(seq1)
  seq2=as.character(seq2)
  if(length(seq1)!=1|length(seq2)!=1){stop("seq1 and seq2 must all contain 1 sequence.");}
  if(nchar(seq1)!=nchar(seq2)|seq1==""){stop("residue numbers in seq1 and seq2 must equal, and can't be 0.");}
  #<<<control the input of seq1 and seq2 end<<<
  #>>>main manipulation begin>>>
  seq1=strsplit(seq1,split="")[[1]] %>% toupper()
  seq2=strsplit(seq2,split="")[[1]] %>% toupper()
  if(unique(c(seq1,seq2)) %>% {!all(. %in% DNA_ALPHABET)}){
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
  ans="A-c";tempo="a-g";
  ans="a-c";tempo="a+g";
  ans="a.c";tempo="agg";
  mergeIUPACLetters.for.2.sequences(seq1=ans,seq2=tempo)
}
