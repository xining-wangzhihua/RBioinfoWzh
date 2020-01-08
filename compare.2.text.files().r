#last edited at 20190108

#a function to export different lines between 2 nearly same documents

warning("\ncompare.2.text.files() can't compare 2 documents with different numbers of lines\n",
        "ask the authour to use setdiff(), match(), etc to revise this function if needed\n\n")

library(Biostrings)
library(R.utils)

compare.2.text.files=function(file1,file2,align_and_trim=TRUE){
  #>>>internal functions begin>>>
  f_trim=function(msu,spbu){
    if(!is.character(msu)|!is.character(spbu)|length(msu)!=1|length(spbu)!=1|msu==spbu){
      stop("invalid input")
    };if(msu==""){msu="-";};if(spbu==""){spbu="-";}
    #
    ans=pairwiseAlignment(pattern=msu,subject=spbu,type="global",gapOpening=10,gapExtension=6)
    msu=strsplit(as.character(Biostrings::pattern(ans)),split="")[[1]]
    spbu=strsplit(as.character(Biostrings::subject(ans)),split="")[[1]]
    rm(ans)
    #
    ans=seqToIntervals( (1:length(msu))[msu!=spbu] )
    if(nrow(ans)==0){stop("internal bug")}
    tempo=matrix(data="",nrow=nrow(ans),ncol=2)
    for(i in 1:nrow(ans)){
      tempo[i,1]=paste0(msu[ans[i,1]:ans[i,2]],collapse="")
      tempo[i,2]=paste0(spbu[ans[i,1]:ans[i,2]],collapse="")
    }
    msu=paste0(tempo[,1],collapse="\t")
    spbu=paste0(tempo[,2],collapse="\t")
    rm(ans,tempo,i)
    #
    return(c(msu,spbu))#rm(msu,spbu)
  }
  #<<<internal functions end<<<
  #>>>control the input of file1 and file2 begin>>>
  if(!is.vector(file1)|is.list(file1)|!is.vector(file2)|is.list(file2)){
    stop("each of file1 and file2 must be a character string.")
  }
  if(length(file1)!=1|length(file2)!=1){
    stop("length(file1) and length(file2) must both be 1.")
  }
  if(!file.exists(file1)|!file.exists(file2)){
    stop("the 2 files specified by file1 and file2 must exist.")
  }
  x=readLines(con=file1);y=readLines(con=file2);
  if(length(x)!=length(y)){stop("numbers of lines in file1 and file2 must equal.")}
  #<<<control the input of file1 and file2 end<<<
  #>>>control the input of align_and_trim begin>>>
  if(length(align_and_trim)!=1){stop("align_and_trim must be one of TRUE or FALSE")}
  align_and_trim=as.logical(align_and_trim)
  if(! align_and_trim %in% c(TRUE,FALSE)){stop("align_and_trim must be one of TRUE or FALSE")}
  #<<<control the input of align_and_trim end<<<
  #>>>main manipulation begin>>>
  tulip=(x!=y)
  if(length(tulip)==0){stop("numbers of lines in file1 and file2 can't be 0")}
  tulip=data.frame(index=(1:length(tulip))[tulip],x=x[tulip],y=y[tulip],stringsAsFactors=FALSE)
  if(nrow(tulip)==0){
    output_file=strftime(x=Sys.time(),format=", %Y%m%d_%H%M%S.txt")
    output_file=paste0("same, ",basename(file1)," & ",basename(file2),output_file,collapse="")
    cat("",file=output_file,append=FALSE)
  }
  if(nrow(tulip)>0 & align_and_trim)for(i in 1:nrow(tulip)){
    tulip[i,2:3]=f_trim(tulip[i,2],tulip[i,3])
  }
  if(nrow(tulip)==1){
    output_file=strftime(x=Sys.time(),format=", %Y%m%d_%H%M%S.txt")
    output_file=paste0("different, ",basename(file1)," & ",basename(file2),output_file,collapse="")
    cat("\nline number: ",tulip$index,"\n",tulip$x,"\n",tulip$y,"\n",file=output_file,append=FALSE)
  }
  if(nrow(tulip)>1){
    output_file=strftime(x=Sys.time(),format=", %Y%m%d_%H%M%S.txt")
    output_file=paste0("different, ",basename(file1)," & ",basename(file2),output_file,collapse="")
    cat("\nline number: ",tulip$index[1],"\n",tulip$x[1],"\n",tulip$y[1],"\n",file=output_file,append=FALSE,sep="")
    for(i in 2:nrow(tulip)){
      cat("\nline number: ",tulip$index[i],"\n",tulip$x[i],"\n",tulip$y[i],"\n",file=output_file,append=TRUE,sep="")
    }
  }
  message("\n\nresult of comparison is written in: ",output_file,"\n\n")
  #<<<main manipulation end<<<
  i=0;rm(file1,file2,tulip,i)
  invisible(output_file)#rm(output_file)
}

if(FALSE){
  #here are the test codes
  file1="~/../OneDrive/Sep13_atrata/New folder/its haplotype with structure.xml"
  file2="~/../OneDrive/Sep13_atrata/New folder/its haplotype without structure.xml"
  compare.2.text.files(file1,file2)
  compare.2.text.files(file1,file2,align_and_trim=FALSE)
  rm(file1,file2,compare.2.text.files)
}
