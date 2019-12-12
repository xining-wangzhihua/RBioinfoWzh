#last edited at 20191212

is.valid.xml.tag.name=function(strings){
  #>>>control the input of strings begin>>>
  if(!is.vector(strings)|is.list(strings)){stop("strings must be a character vector.")}
  strings=as.character(strings)
  #<<<control the input of strings end<<<
  #>>>main manipulation begin>>>
  l=length(strings)
  result=rep(TRUE,times=l)
  result[is.na(strings)]=FALSE
  if(l>0)for(i in 1:l)if(result[i])if(strings[i]==""){
    result[i]=FALSE
  }
  if(l>0)for(i in 1:l)if(result[i])if(grepl(pattern="[^[:alnum:][:punct:]]",x=strings[i])){
    result[i]=FALSE
  }
  if(l>0)for(i in 1:l)if(result[i])if(!grepl(pattern="^[[:alpha:]_]",x=strings[i])){
    result[i]=FALSE
  }
  #<<<main manipulation end<<<
  i=0;rm(strings,l,i)
  return(result)#rm(result)
}

if(FALSE){
  #here are the test codes
  ans=NULL
  ans=list()
  ans=matrix(1:6,nrow=2,ncol=3)
  ans=character()
  ans=integer()
  ans=NA
  ans=c(NA,NULL,NaN)
  ans=c("a","a b","a\tb","a1","1a","a_1","_1a",NA,""," ")
  is.valid.xml.tag.name(strings=ans)
}
