#last edited at 20191212

paste.xml.elements=function(names=character(),attrs=list(),values=list(),pad_start_tag=FALSE){
  x="names, attrs, values must be of same length (or length 1 to be replicated)"
  rm(x)#x serves as annotation
  nms=names;rm(names);
  #>>>r package dependency, function dependency, internal functions begin>>>
  library(package="gtools")
  library(package="stringr")
  library(package="tibble")
  library(package="XML")
  source(file="https://github.com/ywd5/r-zm/raw/master/is.valid.xml.tag.name().r")
  f_entity_reference=function(va){
    #for internal usage only, input validity isn't checked
    va=gsub(pattern="&",replacement="&amp;",x=va,fixed=TRUE)
    va=gsub(pattern="<",replacement="&lt;",x=va,fixed=TRUE)
    va=gsub(pattern=">",replacement="&gt;",x=va,fixed=TRUE)
    va=gsub(pattern="'",replacement="&apos;",x=va,fixed=TRUE)
    va=gsub(pattern="\"",replacement="&quot;",x=va,fixed=TRUE)
    return(va)#rm(va)
  }
  #<<<r package dependency, function dependency, internal functions end<<<
  #>>>control the structure of names, attrs and values begin>>>
  if(!is.vector(nms)|is.list(nms)){stop("names must be a character vector.")}
  if(!is.list(attrs)){stop("attrs must be a list.")}
  if(!is.list(values)){stop("values must be a list.")}
  #
  if(length(nms)==0|length(values)==0){stop("names or values can't be empty.");}
  if(length(attrs)==0){attrs=list(character());}
  l=max(length(nms),length(attrs),length(values))
  if( length(setdiff(x=c(length(nms),length(attrs),length(values)),y=c(1,l)))>0 ){
    stop("length of names, attrs, values must be either 1 or same with others.")
  }
  if(length(nms)==1){nms=rep(nms,times=l);}
  if(length(attrs)==1){attrs=rep(attrs,times=l);}
  if(length(values)==1){values=rep(values,times=l);}
  #<<<control the structure of names, attrs and values end<<<
  #>>>control the content of names, attrs and values begin>>>
  nms=as.character(nms)
  if(!all( is.valid.xml.tag.name(strings=nms) )){stop("some items of names are invalid.")}
  #
  for(i in 1:l)if(length(attrs[[i]])==0){attrs[[i]]=character();}
  for(i in 1:l)if(length(attrs[[i]])!=0){
    ans=attrs[[i]]
    if(!is.vector(ans)|is.list(ans)){stop("attrs[[i]] must be a characters vector.")}
    if(anyNA(ans)){stop("attrs[[i]] can't contain NA.")}
    if(gtools::invalid(names(ans))){stop("names(attrs[[i]]) must be specified.")}
    if(!all( is.valid.xml.tag.name(strings=names(ans)) )){stop("some names(attrs[[i]]) are invalid")}
    attrs[[i]]=setNames(object=as.character(ans),nm=names(ans))
    rm(ans)
  }
  #
  for(i in 1:l)if(length(values[[i]])==0){values[[i]]=character();}
  if( any(lengths(values)>1) ){stop("values[[i]] can't contain more than 1 strings.")}
  if( anyNA(unlist(values)) ){stop("values[[i]] can't be NA.")}
  for(i in 1:l)if(length(values[[i]])==1){
    if(!is.vector(values[[i]])|is.list(values[[i]])){stop("values[[i]] must be a string.")}
    values[[i]]=as.character(values[[i]])
  }
  #<<<control the content of names, attrs and values end<<<
  #>>>control the input of pad_start_tag begin>>>
  pad_start_tag=unname(pad_start_tag)
  if(!identical(pad_start_tag,TRUE) & !identical(pad_start_tag,FALSE)){
    stop("pad_start_tag must be TRUE or FALSE")
  }
  #<<<control the input of pad_start_tag end<<<
  #>>>change 5 entity references in xml begin>>>
  nms=f_entity_reference(va=nms)
  for(i in 1:l)if(length(attrs[[i]])!=0){
    attrs[[i]]=f_entity_reference(va=attrs[[i]])
    names(attrs[[i]])=f_entity_reference(va=names(attrs[[i]]))
  }
  for(i in 1:l)if(length(values[[i]])!=0){
    values[[i]]=f_entity_reference(va=values[[i]])
  }
  #<<<change 5 entity references in xml end<<<
  #>>>main manipulation begin>>>
  result=rep("",times=l)
  if(l>0)for(i in 1:l){
    result[i]=toString( xmlNode(name=nms[i],attrs=attrs[[i]],.children=values[[i]]) )
  }
  if(pad_start_tag){
    result=str_match(string=result,pattern="^(<[^>]+>)(.*</[^>]+>|)$")
    result=as_tibble(result,.name_repair="minimal")[,-1]
    if(anyNA( result[[1]] )){stop("regular expression bug")}
    if(!identical( dim(result),as.integer(c(l,2)) )){stop("bug")}
    result[[1]]=str_pad(string=result[[1]],width=max(nchar(result[[1]])),side="left")
    result=paste0(result[[1]],result[[2]])
  }
  result=paste0(result,"\n")
  #<<<main manipulation end<<<
  rm(nms,attrs,values,pad_start_tag,f_entity_reference,l,i)
  return(result)#rm(result)
}


if(FALSE){
  #here are the test codes
  msu="a";spbu=list("h");itmo=list(c(x="u"));
  msu="a<b";spbu=list("h>i");itmo=list(c(`x&y\'z`="u"));
  msu="a";spbu=list("h");itmo=list(c(x="u"),character());
  msu="a";spbu=list("h",character());itmo=list(c(x="u"));
  msu=c("a","b","ab");spbu=list("h",character(),"hi");itmo=list(c(x="u"),c(`xyz`="uvw"),character());
  msu=c("a","a>b<c&d'e\"f","a-b")
  spbu=list("h","h>i<j&k'l\"m",character())
  itmo=list(c(x="u"),c(x="u",`x>y<z&x'y\"z\\x`="u>v<w&u'v\"w"),character())
  message(paste0( paste.xml.elements(names=msu,attrs=itmo,values=spbu) ,"\n"))
  message(paste0( paste.xml.elements(names=msu,attrs=itmo,values=spbu,pad_start_tag=TRUE) ,"\n"))
  rm(msu,spbu,itmo)
}
