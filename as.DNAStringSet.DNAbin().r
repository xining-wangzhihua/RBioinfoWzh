#a function to convert DNAbin object in ape to DNAStringSet object in Biostrings

#get dependency
if(FALSE){
  library(magrittr)
  library(Biostrings); library(ape);
}

as.DNAStringSet.DNAbin=function(x){
  #control the input of x-------------------------------------------------------
  if(!is(object=x, class2="DNAbin")){stop("x must be a DNAbin object defined in ape.");}
  if(is.matrix(x)){x %<>% as.list();}
  if(!is.list(x)){stop("x must be in list or matrix form, but not vector or other forms");}
  if(length(x) == 0){stop("length(x) can't 0");}
  #check if DNAbin contain invalid characters for DNAStringSet------------------
  lapply(X=x, FUN=as.character) %>% lapply(FUN=anyNA) %>%
    unlist() %>% base::which() %>% paste0(collapse=" ") %>%
    {if(. != ""){stop("the following sequences contain NA: ", .)};}
  lapply(X=x, FUN=as.character) %>% lapply(FUN=function(datum){any(datum == "?")}) %>%
    unlist() %>% base::which() %>% paste0(collapse=" ") %>%
    {if(. != ""){stop("the following sequences contain \"?\": ", .)};}
  #main manipulation------------------------------------------------------------
  no_names <- list(before=is.null(labels(x)), after=FALSE)
  x %<>% base::as.character() %>% lapply(FUN=paste0, collapse="") %>%
    unlist() %>% base::toupper() %>% DNAStringSet()
  no_names$after <- is.null(names(x))
  if(no_names$before != no_names$after){stop("sequence names bug");}
  remove(no_names)
  #return x---------------------------------------------------------------------
  return(x)
}

#here are the test codes
if(FALSE){
  x <- list(c("a","c","g","t","M","r"), "") %>% as.DNAbin()
  x <- list(c("a","c","g","t","M","r"), character()) %>% as.DNAbin()
  x <- list(c("a","c","g","t","M","r"), character()) %>% as.DNAbin() %>% setNames(nm=c("x","y"))
  as.DNAStringSet.DNAbin(x)
}
