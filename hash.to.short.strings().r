#last edited at 20191212

hash.to.short.strings=function(names.arg){
  #>>>r package dependency, function dependency, internal functions begin>>>
  library("tibble")
  #<<<r package dependency, function dependency, internal functions end<<<
  #>>>control the input of names.arg and same_nchar begin>>>
  if(is.null(names.arg)){stop("\"names.arg\" can't be NULL.")}
  if(is.factor(names.arg)){names.arg=as.vector(names.arg);}
  if(!is.vector(names.arg)|is.list(names.arg)){
    stop("\"names.arg\" must be a vector (but not a list).")
  }
  #<<<control the input of names.arg and same_nchar end<<<
  #>>>main manipulation begin>>>
  character_list=c(0:9,letters,LETTERS)
  if(length(names.arg)==0){names.arg=character()}
  if(length(names.arg)!=0)if(length(unique(names.arg))==1){
    names.arg=rep("0",times=length(names.arg))
  }
  if(length(names.arg)!=0)if(length(unique(names.arg))!=1){
    names.arg=match(x=names.arg,table=unique(names.arg))
    l=max(names.arg)
    #
    ans=ceiling( log(x=l,base=length(character_list)) )
    ans=rep(list(character_list),times=ans)
    ans=as_tibble(rev( expand.grid(ans,KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE) ))
    ans=ans[1:l,]
    for(i in 1:l){ans[i,1]=paste0(ans[i,],collapse="");}
    ans=ans[[1]]
    if(anyDuplicated(ans)!=0){stop("bug")}
    #
    names.arg=ans[names.arg]
    rm(l,ans,i)
  }
  #<<<main manipulation end<<<
  return(names.arg)#rm(names.arg,character_list)
}

if(FALSE){
  #here are the test codes
  ans=c("*","(",")","-","_","=","+","\\","|",";",":","/","?")
  msu=rep("",times=9999)
  for(i in 1:length(msu)){msu[i]=paste0(sample(x=ans,size=3,replace=TRUE),collapse="")}
  msu=sample(x=1:999,size=length(msu),replace=TRUE)
  #
  spbu=hash.to.short.strings(names.arg=msu)
  if(!identical( match(msu,unique(msu)),match(spbu,unique(spbu)) )){stop("bug")}
  if(!identical( duplicated(msu),duplicated(spbu) )){stop("bug")}
  anyDuplicated(spbu)
  rm(ans,i,msu,spbu)
}
