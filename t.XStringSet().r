#last edited at 20200101

#a function to transpose an XStringSet object like a matrix.

t.XStringSet=function(xss){
  library(Biostrings)
  #>>>control the input of xss begin>>>
  if(!is(object=xss,class2="XStringSet")){stop("xss must be an XStringSet object.")}
  if(length(xss)==0){stop("length(xss) can't be 0.")}
  if(length(unique(nchar(xss)))!=1){stop("all nchar(xss) must be same.")}
  if(nchar(xss[[1]])==0){stop("sequences in xss can't be empty.")}
  #<<<control the input of xss end<<<
  #>>>main manipulation begin>>>
  xss=unname(xss)
  l=nchar(xss[[1]])
  result=rep(replace(x=xss[1],list=1,values=""),times=l)
  for(i in 1:l){
    result[[i]]=paste0(as.character(subseq(x=xss,start=i,end=i)),collapse="")
  }
  #<<<main manipulation end<<<
  rm(xss,l,i)
  return(result)
}

if(FALSE){
  #here are the test codes
  xss=DNAStringSet(c("acgt","acg"))
  xss=DNAStringSet(c("acgt","acgg"))
  xss=AAStringSet(c("arndc","arndd"))
  xss=BStringSet(c("arndc","arndd"))
  t.XStringSet(xss)
  identical(t.XStringSet(t.XStringSet(xss)),xss)
  all.equal(t.XStringSet(t.XStringSet(xss)),xss,check.attributes=FALSE)
}
