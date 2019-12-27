#last edited at 20191227

#a function to count unique bases occur at each site of aligned sequences

site.polymorphism=function(dnabin){
  #>>>r package dependency, function dependency, internal functions begin>>>
  library(magrittr)
  library(ape)
  #<<<r package dependency, function dependency, internal functions end<<<
  #>>>control the input of dnabin begin>>>
  if(!is(object=dnabin,class2="DNAbin")){stop("dnabin must be a DNAbin object.")}
  dnabin=as.matrix(dnabin)
  if(nrow(dnabin)==0){stop("dnabin must contain at least one sequence.")}
  #<<<control the input of dnabin end<<<
  #>>>main manipulation begin>>>
  result=rep(list(character()),times=ncol(dnabin))
  if(ncol(dnabin)>0)for(i in 1:ncol(dnabin)){
    result[[i]]=as.character(dnabin[,i]) %>% c() %>% unique()
  }
  #<<<main manipulation end<<<
  i=0;rm(dnabin,i)
  return(result)
}

if(FALSE){
  #here are the test codes
  ans=seqinr::s2c("acgtacggacctaagt") %>% rbind(.,.,.)
  ans[1,8]="t";ans[2,11]=NA;ans[3,14]="C";
  ans=as.DNAbin(ans)
  site.polymorphism(dnabin=ans)
  site.polymorphism(dnabin=as.list(ans))
  site.polymorphism(dnabin=ans[,1])
  site.polymorphism(dnabin=ans[FALSE,])
  site.polymorphism(dnabin=ans[,FALSE])
  #a downstream usage:
  ans[,lengths(site.polymorphism(dnabin=ans))!=1] %>% as.character()
}
