#a function to convert AAbin object in ape to AAStringSet object in Biostrings

#get dependency
if(FALSE){
  library(magrittr)
  library(Biostrings); library(ape);
}

as.AAStringSet.AAbin <- function(x){
  #control the input of x-------------------------------------------------------
  if(!is(object=x, class2="AAbin")){stop("x must be an AAbin object in ape");}
  if(is.matrix(x)){x %<>% as.list();}
  if(!is.list(x)){stop("x must be in list or matrix form, but not vector or other forms");}
  if(length(x) == 0){stop("length(x) can't 0");}
  #check if AAbin contain invalid characters for AAStringSet--------------------
  lapply(X=x, FUN=as.character) %>% lapply(FUN=anyNA) %>%
    unlist() %>% base::which() %>% paste0(collapse=" ") %>%
    {if(. != ""){stop("the following sequences contain NA: ", .)};}
  lapply(X=x, FUN=as.character) %>% lapply(FUN=function(datum){any(datum == "?")}) %>%
    unlist() %>% base::which() %>% paste0(collapse=" ") %>%
    {if(. != ""){stop("the following sequences contain \"?\": ", .)};}
  #main manipulation------------------------------------------------------------
  no_names <- list(before=is.null(labels(x)), after=FALSE)
  x %<>% base::as.character() %>% lapply(FUN=paste0, collapse="") %>%
    unlist() %>% base::toupper() %>% AAStringSet()
  no_names$after <- is.null(names(x))
  if(no_names$before != no_names$after){stop("sequence names bug");}
  remove(no_names)
  #return x---------------------------------------------------------------------
  return(x)
}

#here are the test codes
if(FALSE){
  #case 1
  x <- list(sample(AA_ALPHABET, size=8), sample(AA_ALPHABET, size=9), character(), "")
  names(x) <- c("a1", "a2", "a3", "a4")[1:length(x)]
  #case 2
  x <- list(sample(AA_ALPHABET, size=9), sample(AA_ALPHABET, size=9))
  names(x) <- c("a1", "a2", "a3", "a4")[1:length(x)]
  #case 3
  x <- list()
  #case 4
  x <- list(a1=c("X", "z"), a2=c("g", "?", "Q"))
  #
  as.AAbin(x) %>% as.character()
  as.AAbin(x) %>% as.AAStringSet.AAbin()
  as.AAbin(x) %>% .[[1]] %>% as.AAStringSet.AAbin()
  as.AAbin(x) %>% as.matrix() %>% as.AAStringSet.AAbin()
}
